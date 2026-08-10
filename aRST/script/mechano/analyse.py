from aRST.script.exploration.intraD import dissoxtbscan,dissombff
import aRST.script.initialization as storage
from aRST.script.exploration.structure import StrucGraph
from aRST.script.do_scan import get_endxyz_from_mb
from aRST.script.toolbox import create_folder,del_dir
from aRST.script.exploration.pathsearch import SinglePathSearch
from aRST.script.head import print_simulation_details, print_jobhead
from aRST.script.exploration.reactionsearch import PathEval
from aRST.script.exploration.network import rxnnetwork

from os.path import join
import networkx as nx
import numpy as np
import shutil
import json
import copy


# class MechanophoreTreatment:
#     def __init__(self,structure)


def get_streched_geom(struc, entry, count=0):
    output_d = {}

    scanwd = storage.setting.searching_setting.general.scanwd
    s_general = storage.allstruc.get(struc).general
    charge = s_general.charge
    nel = s_general.nel
    uhf_l = [0, 2] if nel % 2 == 0 else [1]
    disso_setting = storage.setting.searching_setting.disso
    stretch_scale = disso_setting.stretch_scale

    coordinate = s_general.coordinate
    ele_l = s_general.ele_l
    # get atomic connection
    atom_connection = copy.deepcopy(s_general.atomic_connection)
    atom_1 = entry[0]
    atom_2 = entry[1]
    atom_1_, atom_2_ = atom_1 - 1, atom_2 - 1
    at_l = [[atom_1_, atom_2_]]
    mbo = entry[2]
    r0 = entry[3]

    target_dist = r0 * stretch_scale
    atom_connection.append({"atoms": [atom_1, atom_2], "value": target_dist})

    mbxyz0, ele_l0 = get_endxyz_from_mb(atom_connection, coordinate, ele_l, set_angles=False)
    mbxyz, mbele_l = get_endxyz_from_mb(atom_connection, mbxyz0, ele_l0, set_angles=True)

    # build up push_l for calling SinglePathSearch
    push_l = [[atom_1, atom_2, r0, target_dist]]

    rc = 1
    for u in uhf_l:
        scan_wd = join(scanwd, f"0_{count}","rupture_force")
        create_folder(scan_wd)
        reactantinfo = {"reactants": [struc],
                        "unreacted": [],
                        "at_l": at_l,
                        "simuhf": u,
                        "mbo": mbo,
                        "scanwd": f"0_{count}",
                        "coordinate": coordinate,
                        "ele_l": ele_l,
                        "mb": s_general.mb
                        }  # reactants info
        results, status = SinglePathSearch(mbxyz,
                                           mbele_l,
                                           u, charge, scan_wd,
                                           reactants=[struc]).run_xtbopt_for_mbff_search(
            push_l=push_l, rxn_type="dissoD", reactantinfo=reactantinfo, nel=nel)
        if status == 0:
            rc = 0
            # move wbo and charge to xtb folder for later read for single-struc process
            tmp = create_folder(join(scanwd, f"0_{count}","rupture_force", "xtb"))
            shutil.copyfile(join(scanwd, f"0_{count}","rupture_force", "wbo"), join(tmp,"wbo"))
            shutil.copyfile(join(scanwd, f"0_{count}","rupture_force", "charges"), join(tmp,"charges"))
        if results:
            output_d[f"0_{count}"] = results

        count += 1
    return output_d, count, rc


def follow_up_process(entry):
    # check number of fragments
    product_info = entry[1]
    pcoordinate = product_info["coordinate"]
    pele_l = product_info["ele_l"]
    ptopo_l = product_info["atomic_connection"]
    puhf = product_info["uhf"]
    pcharge = product_info["charge"]

    reactant_info = entry[0]
    scanwd = reactant_info["scanwd"]
    at_l = reactant_info["at_l"]

    nfragments = nx.number_connected_components(StrucGraph(pcoordinate, pele_l, ptopo_l).graph)
    if nfragments > 1:
        print(f"Warning: Rupture of bond {at_l}, uhf={puhf} lead to multipule fragments!")

    # now get the E_sp of constrained geom
    comdir = join(storage.setting.searching_setting.general.scanwd, scanwd,"rupture_force")
    create_folder(comdir)
    features = storage.Singlestruc(pcoordinate, pele_l, pcharge, puhf + 1, ).process_workflow(comdir,
                                                                                              ideal=False,
                                                                                              opt=False)
    processed_entry = list(entry)
    processed_entry.append(features)
    return processed_entry


def FmaxAnalyse(reactants_ini, ff, debugtrigger=False):
    if len(reactants_ini) > 1:
        print(
            "Warning: More than one initial structure was provided. Mechanical analysis will be performed only for the first molecule.")
    struc = reactants_ini[0]
    print(f"Start doing mechanical analysis with given structues: ", struc)
    print()

    s_general = storage.allstruc.get(struc).general
    coordinate = s_general.coordinate
    ele_l = s_general.ele_l
    topo_l = s_general.atomic_connection
    nel = s_general.nel

    # Estimate the potential rupture sites
    ## check if ext define topo search range
    general_setting = storage.setting.searching_setting.general
    disso_setting = storage.setting.searching_setting.disso
    terminal_pair = disso_setting.terminal_pair
    ncandidates = disso_setting.number
    if terminal_pair:
        print(f"Force application nodes have been defined as {terminal_pair}")
        terminal_pair_1 = [(a + 1, b + 1) for [a, b] in terminal_pair]
        obj = StrucGraph(coordinate, ele_l, topo_l, ftp_analyse=True, terminal_pair=terminal_pair_1)
    else:
        obj = StrucGraph(coordinate, ele_l, topo_l, ftp_analyse=True, terminal_pair=None)

    print("Following force tranmission paths are considered:")
    for bb in obj.all_backbones:
        print([_ - 1 for _ in bb])
    print()

    # Estimate the potential rupture sites
    rupture_bonds_all = obj.get_assess_bond_set()
    # write
    with open(join(general_setting.wd, "rupture_sites_ranking"), "w") as f:
        for item in rupture_bonds_all:
            f.write(json.dumps(item) + "\n")

    rupture_bonds = rupture_bonds_all[:ncandidates]

    # get_streched geom
    candidate_results = {}
    searchcount = 0
    for entry in rupture_bonds:
        res, searchcount, rc = get_streched_geom(struc, entry, count=searchcount)
        candidate_results.update(res)

    # do follow-up prosee for streched geom
    processed_candidate_results = {}
    for scan_name, result in candidate_results.items():
        processed_result = follow_up_process(result)
        processed_candidate_results[scan_name] = processed_result

    # for the same at_l, only keep the streched res with lower E_sp
    results_filtered = {}
    for scan_name, result in processed_candidate_results.items():
        at_l_hash = tuple(tuple(x) for x in result[0]["at_l"])
        esp = result[3]["Esp_gas"]
        if (
                at_l_hash not in results_filtered
                or esp < results_filtered[at_l_hash][3]["Esp_gas"]
        ):
            results_filtered[at_l_hash] = result

    # get Fmax
    out = []
    s_energy = storage.allstruc.get(struc).energy
    E_0 = s_energy.Esp_gas

    if disso_setting.use_d0:
        E_0 += s_energy.ZPE
    for result in results_filtered.values():
        R_info = result[0]
        at_l = R_info["at_l"][0]
        target = {at_l[0], at_l[1]}
        mask = ff["atoms"].apply(lambda x: set(map(int, str(x).strip("[]").split())) == target)
        if not mask.any():
            print(f"Warning: given FF file doesn't contains enrty for atl={at_l}")
            continue

        Pstreched_info = result[3]
        E_1 = Pstreched_info["Esp_gas"]
        if disso_setting.use_d0:
            E_1 += Pstreched_info["ZPE"]
        delta_E_eh = E_1 - E_0
        parameter = ff.loc[mask, "parameter"].iloc[0]
        k = 0.5 * (parameter ** 2)
        r0 = ff.loc[mask, "reference_value"].iloc[0]
        a = np.sqrt(k / (2 * delta_E_eh))  # unit: 1/bohr
        Fmax = a * delta_E_eh / 2
        Fmax_nN = Fmax * 82.387  # EH_PER_BOHR_TO_NN = 82.387

        mbo = R_info["mbo"]
        uhf = R_info["simuhf"]
        row = {
            "atl": at_l,
            "mbo": mbo,
            "uhf": uhf,
            "deltaE[eh]": delta_E_eh,
            "Fmax[nN]": Fmax_nN,
            "k": np.round(k, 5),
            "a[1/bohr]": np.round(a, 6),
        }
        out.append((result, row))
    return out


def RuptureProductsAnalyse(out,debugtrigger=False):
    searchcount = 0
    scanmode = storage.setting.searching_setting.general.scan_mode
    forcescanwd = "rupture_product"
    output_D = {}
    for item in out:
        results = item[0]
        row = item[1]
        # make entry for explore func
        entry = {"reactants": [0], "atl": [row['atl']], "mbo": row['mbo']}

        if scanmode==0:
            res, searchcount, rc = dissoxtbscan([0], entry, "0", searchcount,forcescanwd=forcescanwd,checkpoint=False)

        else:
            res, searchcount, rc = dissombff([0], entry, "0", searchcount, forcescanwd=forcescanwd,checkpoint=False)
        output_D.update(res)

    print_info = []
    for item in output_D.values():
        reactants_info = item[0]
        endxyz_info = item[1]
        reaction_info = item[2]
        rxninfo = PathEval(reactants_info, endxyz_info, reaction_info, debugtrigger).run()
        print_info.append(rxninfo)

    tseval = storage.setting.searching_setting.general.tseval
    evalmode = storage.setting.searching_setting.general.evalmode
    print_info = rxnnetwork.self_check(print_info)
    print_simulation_details(print_info, "0", tseval=tseval, evalmode=evalmode)
    del_dir(not debugtrigger, storage.setting.searching_setting.general.bufferwd)
    return print_info
