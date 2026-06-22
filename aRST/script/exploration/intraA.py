# Reaction simulation for intramolecular association
import json
from os.path import join
import numpy as np
from aRST.script.spr import SPR
from aRST.script.do_scan import get_endxyz_from_mb, build_push_l
from aRST.script.exploration.pathsearch import SinglePathSearch
import aRST.script.initialization as storage
from aRST.script.head import print_react_l
from aRST.script.toolbox import create_folder,write_xyz
import copy


def get_intraA_react_l(reactants0, writefile=True):
    generalsetting = storage.setting.searching_setting.general
    intrasetting = storage.setting.searching_setting.intra
    react_l = SPR(reactants0,
                  generalsetting.spr_mode,
                  ).intra(
        equivalent_afogap_tol=generalsetting.equivalent_afogap_tol,
        exclude_elements=generalsetting.exclude_elements,
        exclude_pairs=generalsetting.exclude_pairs if generalsetting.cycle == 1 else None
    )
    if writefile:
        print_react_l(reactants0, react_l)
        # write full react_l in file
        with open(join(generalsetting.wd, "react_l"), "w") as f:
            for item in react_l:
                f.write(json.dumps(item) + "\n")

    cut = 2 if intrasetting.concerted_pair else 1
    for item in react_l:
        pairs = item["atl"][:cut]
        item["atl"] = pairs



    atom_pairs = intrasetting.atom_pair
    if atom_pairs:
        # overwrite react_l with user's input
        # predifined reaction, used in one cycle
        at1, at2 = atom_pairs[0]
        react_l = [{"reactants": [reactants0[0]],
                    "atl": [[at1, at2]],
                    "gap": 0,
                    "force_combination": True}]
        if len(atom_pairs) > 1:
            intrasetting.concerted_pair = True
            for atom_pair in atom_pairs[1:]:
                atl_list = react_l[0]["atl"]
                at3, at4 = atom_pair
                atl_list.append([at3, at4])

    return react_l


def check_scan_uhf_intraA(entry):
    struc = entry["reactants"][0]
    nel = storage.allstruc.get(struc).general.nel
    if nel % 2 == 0:
        return [0, 2]
    elif nel % 2 != 0:
        return [1]


def intraxtbscan(reactants0, entry, reactionpath, count=0):
    generalsetting = storage.setting.searching_setting.general
    xtbsetting = storage.setting.xtb_setting
    scanwd = generalsetting.scanwd
    output_d = {}
    uhf_l = check_scan_uhf_intraA(entry)
    struc = entry["reactants"][0]
    s = storage.allstruc.get(struc)
    coordinate = s.general.coordinate
    ele_l = s.general.ele_l
    at_l = entry["atl"]
    push_l = build_push_l(coordinate, ele_l, at_l)
    charge = s.general.charge
    nel = s.general.nel
    m = s.general.multiplicity

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

    unreacted = sorted((x for x in reactants0 if x != struc), key=int, reverse=True)

    rc = 1
    for u in uhf_l:
        # print(f"Simulating intraA path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...", end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")
        reactantinfo = {"reactants": [struc],
                        "unreacted": unreacted,
                        "at_l": at_l, "simuhf": u,
                        "scanwd": f"{reactionpath}_{count}",
                        "coordinate": coordinate,
                        "ele_l": ele_l,
                        "mb": s.general.mb
                        }  # reactants info
        results, status = SinglePathSearch(coordinate,
                                           ele_l,
                                           u, charge, scan_wd,
                                           reactants=[struc]).run_one_scan_search(constrain_setting,
                                                                                  "intraA", reactantinfo, nel)
        if status == 0:
            rc = 0
        if results:
            output_d[f"{reactionpath}_{count}"] = results
        #     print(f"done.")
        # else:
        #     print(f"failed.")
        count += 1
    return output_d, count, rc


def intrambff(reactants0, entry, reactionpath, count=0):
    generalsetting = storage.setting.searching_setting.general
    xtbsetting = storage.setting.xtb_setting
    scanwd = generalsetting.scanwd
    output_d = {}
    uhf_l = check_scan_uhf_intraA(entry)
    struc = entry["reactants"][0]
    s = storage.allstruc.get(struc)
    coordinate = s.general.coordinate
    ele_l = s.general.ele_l
    at_l = entry["atl"]
    push_l = build_push_l(coordinate, ele_l, at_l)
    charge = s.general.charge
    nel = s.general.nel
    m = s.general.multiplicity

    # get atomic connection (COPY!)
    atom_connection = copy.deepcopy(s.general.atomic_connection)
    for item in push_l:
        at1, at2 = int(item[0]), int(item[1])
        value = np.round(float(item[3]), 3)
        atom_connection.append({"atoms": [at1, at2], "value": value})

    unreacted = sorted((x for x in reactants0 if x != struc), key=int, reverse=True)
    mbxyz0, ele_l0 = get_endxyz_from_mb(atom_connection, coordinate, ele_l, set_angles=False)
    mbxyz, mbele_l = get_endxyz_from_mb(atom_connection, mbxyz0, ele_l0, set_angles=True)

    rc = 1
    for u in uhf_l:
        # print(f"Simulating intraA path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...", end="")
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
                                           reactants=[struc]).run_xtbopt_for_mbff_search(
            push_l=push_l, rxn_type="intraA", reactantinfo=reactantinfo, nel=nel)
        if status == 0:
            rc = 0
        if results:
            output_d[f"{reactionpath}_{count}"] = results
            # print(f"done.")

        # else:
        #     print(f"failed.")
        count += 1
    return output_d, count, rc
