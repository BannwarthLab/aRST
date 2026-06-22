# Reaction simulation for intramolecular association
import json
from os.path import join
import numpy as np
from aRST.script.spr import SPR
from aRST.script.geom import alignment
from aRST.script.do_scan import get_endxyz_from_mb, build_push_l
from aRST.script.exploration.pathsearch import SinglePathSearch
import aRST.script.initialization as storage
from aRST.script.head import print_react_l
from aRST.script.toolbox import create_folder,write_xyz,get_topo_covaradii
from aRST.script.callsys import CallMolBar
from aRST.script.exploration.structure import check_target_bond_change
import copy


def get_interA_react_l(reactants0, writefile=True):
    generalsetting = storage.setting.searching_setting.general
    intersetting = storage.setting.searching_setting.inter
    react_l = SPR(reactants0, generalsetting.spr_mode).inter(
        equivalent_afogap_tol=generalsetting.equivalent_afogap_tol,
        exclude_elements=generalsetting.exclude_elements,
        exclude_pairs=generalsetting.exclude_pairs if generalsetting.cycle == 1 else None
    )
    if writefile:
        # print top 10 react_l
        print_react_l(reactants0, react_l)
        # write full react_l in file
        with open(join(generalsetting.wd, "react_l"), "w") as f:
            for item in react_l:
                f.write(json.dumps(item) + "\n")

    cut = 2 if intersetting.concerted_pair else 1
    for item in react_l:
        pairs = item["atl"][:cut]
        item["atl"] = pairs


    atom_pairs = intersetting.atom_pair
    if atom_pairs:
        # overwrite react_l with user's input
        # predifined reaction, used in one cycle
        at1, at2 = atom_pairs[0]
        react_l = [{"reactants": [reactants0[0], reactants0[1]],
                    "atl": [[at1, at2]],
                    "gap": 0,
                    "force_combination": True}]
        if len(atom_pairs) > 1:
            intersetting.concerted_pair = True
            for atom_pair in atom_pairs[1:]:
                at_combi_list = react_l[0]["atl"]
                at3, at4 = atom_pair
                at_combi_list.append([at3, at4])
    return react_l


def generate_reactive_complex(entry):
    intersetting = storage.setting.searching_setting.inter
    # prepare scan_xyz
    struc1, struc2 = entry["reactants"]
    s1 = storage.allstruc.get(struc1)
    coordinate1 = s1.general.coordinate
    ele_l1 = s1.general.ele_l
    nelA = s1.general.nel

    s2 = storage.allstruc.get(struc2)
    coordinate2 = s2.general.coordinate
    ele_l2 = s2.general.ele_l
    nelB = s2.general.nel
    nat_coord1 = len(ele_l1)

    # generate reactiv complex
    at_l = entry["atl"]
    coordinate_merge, ele_merged = alignment(coordinate1, ele_l1, coordinate2, ele_l2, at_l,
                                             extra_setting={"maxovlap": intersetting.maxovlap,
                                                            "aligning_distance": intersetting.aligning_distance})
    # prepare scan setteing
    charge = sum([k.general.charge for k in [s1, s2]])

    force_combination = False if "force_combination" not in entry.keys() else entry["force_combination"]
    push_l = build_push_l(coordinate_merge, ele_merged, at_l, nat_coord1, force_combination=force_combination
                          )
    return struc1, struc2, at_l, coordinate_merge, ele_merged, nat_coord1, charge, nelA + nelB, push_l


def check_scan_uhf_interA(entry):
    struc1, struc2 = entry["reactants"]
    nelA = storage.allstruc.get(struc1).general.nel
    nelB = storage.allstruc.get(struc2).general.nel
    if nelA % 2 == 0 and nelB % 2 == 0:
        return [0]
    elif (nelA + nelB) % 2 != 0:
        return [1]
    else:
        return [0, 2]


def interxtbscan(reactants0, entry, reactionpath, count=0):
    generalsetting = storage.setting.searching_setting.general
    intersetting = storage.setting.searching_setting.inter
    xtbsetting = storage.setting.xtb_setting
    scanwd = generalsetting.scanwd
    output_d = {}
    # preset output
    uhf_l = check_scan_uhf_interA(entry)
    struc1, struc2, at_l, \
    coordinate_merge, ele_merged, nat_coord1, \
    charge, totnel, push_l = generate_reactive_complex(entry)

    # preparing scan inp
    constrain_setting = {
        "constrain": {
            "force_constant": xtbsetting.k
        },
        "scan": {
            "mode": "concerted" if len(at_l) == 2 else None,
            "distance": push_l if push_l else None,
            "step": xtbsetting.step
        }
    }
    unreacted = sorted((x for x in reactants0 if x != struc1 and x != struc2), key=int, reverse=True)
    aligned_mb = CallMolBar(coordinate_merge, ele_merged).get_mb()

    rc = 1
    for u in uhf_l:
        # print(f"Simulating interA path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...", end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")

        reactantinfo = {"reactants": [struc1, struc2],
                        "unreacted": unreacted,
                        "at_l": at_l, "simuhf": u,
                        "scanwd": f"{reactionpath}_{count}",
                        "coordinate": coordinate_merge,
                        "ele_l": ele_merged,
                        "mb": aligned_mb
                        }  # reactants info
        results, status = SinglePathSearch(coordinate_merge,
                                           ele_merged,
                                           u, charge, scan_wd,
                                           nat_coord1, [struc1, struc2]).run_one_scan_search(
            constrain_setting, "interA", reactantinfo, totnel
        )
        if status == 0:
            rc = 0
        if results:
            output_d[f"{reactionpath}_{count}"] = results

        # else:
        #     print(f"failed.")
        count += 1
    return output_d, count, rc


def intermbff(reactants0, entry, reactionpath, count=0):
    generalsetting = storage.setting.searching_setting.general
    intersetting = storage.setting.searching_setting.inter
    scanwd = generalsetting.scanwd
    output_d = {}
    uhf_l = check_scan_uhf_interA(entry)

    struc1, struc2, at_l, \
    coordinate_merge, ele_merged, nat_coord1, \
    charge, totnel, push_l = generate_reactive_complex(entry)

    # get atomic connection
    atom_connection1 = copy.deepcopy(storage.allstruc.get(struc1).general.atomic_connection)
    atom_connection2 = copy.deepcopy(storage.allstruc.get(struc2).general.atomic_connection)
    add_atom_connection2 = atom_connection2
    for item in add_atom_connection2:
        for at in range(len(item["atoms"])):
            item["atoms"][at] += nat_coord1

    atom_connection = atom_connection1 + add_atom_connection2
    aligned_mb =  CallMolBar(coordinate_merge, ele_merged).get_mb()
    for item in push_l:
        at1, at2 = int(item[0]), int(item[1])
        value = np.round(float(item[3]), 2)
        atom_connection.append({"atoms": [at1, at2], "value": value})

    unreacted = sorted((x for x in reactants0 if x != struc1 and x != struc2), key=int, reverse=True)
    # # test mb set_angles=True can connecting fragments
    # mbxyz, mbele_l = get_endxyz_from_mb(atom_connection, coordinate_merge, ele_merged, set_angles=True)
    # tmp_topo_l =  get_topo_covaradii(mbxyz, mbele_l)
    # maped_atl = [[item[0], item[1]] for item in push_l]
    # target_change = check_target_bond_change(maped_atl, tmp_topo_l, "interA")
    # if not target_change:
    mbxyz0, ele_l0 = get_endxyz_from_mb(atom_connection, coordinate_merge, ele_merged, set_angles=False)
    atom_connection = get_topo_covaradii(mbxyz0,ele_l0)
    mbxyz, mbele_l = get_endxyz_from_mb(atom_connection, mbxyz0, ele_l0, set_angles=True)
    rc = 1
    for u in uhf_l:
        # print(f"Simulating interA path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...", end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")
        create_folder(scan_wd)
        write_xyz(coordinate_merge, ele_merged, "coord")
        # write_xyz(mbxyz0, ele_l0, "a")
        reactantinfo = {"reactants": [struc1, struc2],
                        "unreacted": unreacted,
                        "at_l": at_l, "simuhf": u,
                        "scanwd": f"{reactionpath}_{count}",
                        "coordinate": coordinate_merge,
                        "ele_l": ele_merged,
                        "mb":aligned_mb
                        }  # reactants info
        results, status = SinglePathSearch(mbxyz,
                                           mbele_l,
                                           u, charge, scan_wd,
                                           nat_coord1,
                                           [struc1, struc2]).run_xtbopt_for_mbff_search(
            push_l=push_l,rxn_type="interA",reactantinfo=reactantinfo,nel=totnel)

        if status == 0:
            rc = 0
        if results:
            output_d[f"{reactionpath}_{count}"] = results

        #     print(f"done.")
        #
        # else:
        #     print(f"failed.")
        count += 1
    return output_d, count, rc
