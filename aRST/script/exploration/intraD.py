import json
from os.path import join
import numpy as np
from aRST.script.spr import SDP
from aRST.script.do_scan import get_endxyz_from_mb
from aRST.script.exploration.pathsearch import SinglePathSearch
import aRST.script.initialization as storage
from aRST.script.head import print_react_l
from aRST.script.geom import Coord
from aRST.script.exploration.structure import StrucGraph
from aRST.script.toolbox import create_folder,write_xyz
import copy


def get_intraD_react_l(reactants0, gaplim=5, writefile=True):
    general_setting = storage.setting.searching_setting.general
    disso_setting = storage.setting.searching_setting.disso
    react_l = SDP(reactants0,
                  exclude_elements=general_setting.exclude_elements,
                  exclude_pairs=general_setting.exclude_pairs if general_setting.cycle == 1 else None
                  ).build_react_l(tol=1e-5)
    cut = 2 if disso_setting.concerted_pair else 1
    for item in react_l:
        pairs = item["atl"][:cut]
        item["atl"] = pairs

    if writefile:
        print_react_l(reactants0, react_l)
        # write full react_l in file
        with open(join(general_setting.wd, "react_l"), "w") as f:
            for item in react_l:
                f.write(json.dumps(item) + "\n")
    atom_pairs = disso_setting.atom_pair
    if atom_pairs:
        # overwrite react_l with user's input
        # predifined reaction, used in one cycle
        at1, at2 = atom_pairs[0]
        react_l = [{"reactants": [reactants0[0]],
                    "atl": [[at1, at2]],
                    "gap": 0,
                    "force_combination": True}]
        if len(atom_pairs) > 1:
            disso_setting.concerted_pair = True
            for atom_pair in atom_pairs[1:]:
                at_combi_list = react_l[0]["atl"]
                at3, at4 = atom_pair
                at_combi_list.append([at3, at4])
    return react_l


def check_scan_uhf_D(entry):
    struc = entry["reactants"][0]
    nel = storage.allstruc.get(struc).general.nel
    if nel % 2 == 0:
        return [0, 2]
    elif nel % 2 != 0:
        return [1]


def determine_pulling_distance(entry):
    disso_setting = storage.setting.searching_setting.disso
    struc = entry["reactants"][0]
    s_general = storage.allstruc.get(struc).general
    coordinate = s_general.coordinate
    ele_l = s_general.ele_l
    at_l = entry["atl"]

    instance = Coord(coordinate, ele_l)

    # output
    push_l = []

    # now determine the definition of pulling distance
    if disso_setting.pulling_distance_add:
        ## A: Defined as an adding distance to the current distance
        add_dis = disso_setting.pulling_distance_add
        for pair in at_l[:2]:
            at1, at2 = map(int, pair)
            dis_current = np.round(instance.get_ats_distance(at1, at2), 3)
            dis_tar = dis_current + add_dis
            at1_reindex = at1 + 1
            at2_reindex = at2 + 1
            push_l.append([at1_reindex, at2_reindex, dis_current, dis_tar])

    else:
        mul_dis = disso_setting.pulling_distance

        if mul_dis == "auto" and len(at_l) == 1:
            strucgraph = StrucGraph(coordinate, ele_l, s_general.atomic_connection)

            at1, at2 = map(int, at_l[0])  # 0-based
            at1_, at2_ = at1 + 1, at2 + 1  # 1-based for graph

            dis_current = np.round(instance.get_ats_distance(at1, at2), 3)

            if not strucgraph.bond_exist([at1_, at2_]):
                atomic_path, atomic_length = strucgraph.get_atomic_path_of_2atoms(at1_, at2_)

                if atomic_length is not None:
                    push_l.append([at1_, at2_, dis_current, atomic_length])
                    return push_l

            # fallback if bonded OR no atomic path
            mul_dis = 5
        else:
            mul_dis = 5

        for pair in at_l[:2]:
            at1, at2 = map(int, pair)

            dis_current = np.round(instance.get_ats_distance(at1, at2), 3)
            cradii_sum = float(instance.rcov_l[at1] + instance.rcov_l[at2])

            at1_reindex = at1 + 1
            at2_reindex = at2 + 1

            push_l.append([
                at1_reindex,
                at2_reindex,
                dis_current,
                cradii_sum * float(mul_dis)
            ])

        return push_l


def dissoxtbscan(reactants0, entry, reactionpath, count=0):
    scanwd = storage.setting.searching_setting.general.scanwd
    output_d = {}
    uhf_l = check_scan_uhf_D(entry)

    struc = entry["reactants"][0]
    s_general = storage.allstruc.get(struc).general
    coordinate = s_general.coordinate
    ele_l = s_general.ele_l
    at_l = entry["atl"]

    # prepare scan
    charge = s_general.charge
    nel = s_general.nel
    m = s_general.multiplicity

    # preparing scan inp
    push_l = determine_pulling_distance(entry)
    xtbsetting = storage.setting.xtb_setting
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
    unreacted = sorted((x for x in reactants0 if x != struc), key=int, reverse=True)
    rc = 1
    for u in uhf_l:
        # print(f"Simulating disso path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...", end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")
        reactantinfo = {"reactants": [struc],
                        "unreacted": unreacted,
                        "at_l": at_l, "simuhf": u,
                        "scanwd": f"{reactionpath}_{count}",
                        "coordinate": coordinate,
                        "ele_l": ele_l
                        }  # reactants info

        results, status = SinglePathSearch(coordinate,
                                           ele_l,
                                           u, charge, scan_wd,
                                           reactants=[struc]).run_one_scan_search(constrain_setting, "dissoD",
                                                                                  reactantinfo,
                                                                                  nel)
        if status == 0:
            rc = 0
        if results:
            output_d[f"{reactionpath}_{count}"] = results

        # else:
        #     print(f"failed.")
        count += 1
    return output_d, count, rc


def dissombff(reactants0, entry, reactionpath, count=0):
    scanwd = storage.setting.searching_setting.general.scanwd
    output_d = {}
    uhf_l = check_scan_uhf_D(entry)

    struc = entry["reactants"][0]
    s_general = storage.allstruc.get(struc).general
    coordinate = s_general.coordinate
    ele_l = s_general.ele_l
    at_l = entry["atl"]

    # prepare scan
    charge = s_general.charge
    nel = s_general.nel
    m = s_general.multiplicity

    # preparing scan inp
    push_l = determine_pulling_distance(entry)

    # get atomic connection
    atom_connection = copy.deepcopy(s_general.atomic_connection)
    for item in push_l:
        at1 = int(item[0])
        at2 = int(item[1])
        value = np.round(float(item[3]), 2)
        atom_connection.append({"atoms": [at1, at2], "value": value})
    unreacted = sorted((x for x in reactants0 if x != struc), key=int, reverse=True)

    mbxyz0, ele_l0 = get_endxyz_from_mb(atom_connection, coordinate, ele_l, set_angles=False)
    mbxyz, mbele_l = get_endxyz_from_mb(atom_connection, mbxyz0, ele_l0, set_angles=True)

    rc = 1
    for u in uhf_l:
        # print(f"Simulating disso path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...", end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")
        create_folder(scan_wd)
        write_xyz(coordinate, ele_l, "coord")
        reactantinfo = {"reactants": [struc],
                        "unreacted": unreacted,
                        "at_l": at_l, "simuhf": u,
                        "scanwd": f"{reactionpath}_{count}",
                        "coordinate": coordinate,
                        "ele_l": ele_l,
                        "mb":s.general.mb
                        }  # reactants info
        results, status = SinglePathSearch(mbxyz,
                                           mbele_l,
                                           u, charge, scan_wd,
                                           reactants=[struc]).run_xtbopt_for_mbff_search(atom_connection, push_l,
                                                                                         "dissoD", reactantinfo, nel)

        if status == 0:
            rc = 0
        if results:
            output_d[f"{reactionpath}_{count}"] = results

            # print(f"done.")

        # else:
        # print(f"failed.")
        count += 1
    return output_d, count, rc
