# Reaction simulation for intramolecular association
import json
from os.path import join
import shutil
import numpy as np
import pandas as pd
from aRST.script.toolbox import create_folder, write_xyz, read_coord, read_bo, read_atcharge, get_topo
from aRST.script.do_scan import build_push_l, get_endxyz_from_mb
from aRST.script.spr import SPR, SDP
from aRST.script.geom import Coord
from aRST.script.callsys import CallxTB
import aRST.script.initialization as storage
from aRST.script.head import print_react_l
import copy


def get_intraA_react_l(reactants0):
    react_l = SPR(reactants0,
                  storage.setting.searching_setting.general.spr_mode,
                  ).intra(
        equivalent_afogap_tol=storage.setting.searching_setting.general.equivalent_afogap_tol,
        exclude_elements=storage.setting.searching_setting.general.exclude_elements,
        exclude_pairs=storage.setting.searching_setting.general.exclude_pairs
    )
    cut = 2 if storage.setting.searching_setting.intra.concerted_pair else 1
    for item in react_l:
        pairs = item["at_combi"][:cut]
        item["at_combi"] = pairs

    print_react_l(reactants0, react_l)
    # write full react_l in file
    with open(join(storage.setting.searching_setting.general.wd, "react_l"), "w") as f:
        for item in react_l:
            f.write(json.dumps(item) + "\n")

    atom_pairs = storage.setting.searching_setting.intra.atom_pair
    if atom_pairs:
        # overwrite react_l with user's input
        # predifined reaction, used in one cycle
        at1, at2 = atom_pairs[0]
        react_l = [{"struc_combi": [reactants0[0], reactants0[1]],
                    "at_combi": [[at1, at2]],
                    "gap": 0,
                    "force_combination": True}]
        if len(atom_pairs) > 1:
            storage.setting.searching_setting.intra.concerted_pair = True
            for atom_pair in atom_pairs[1:]:
                at_combi_list = react_l[0]["at_combi"]
                at3, at4 = atom_pair
                at_combi_list.append([at3, at4])

    return react_l


def check_scan_uhf_intraA(entry):
    struc = entry["struc"]
    nel = storage.allstruc.get(struc).struc_info.general.nel
    if nel % 2 == 0:
        return [0, 2]
    elif nel % 2 != 0:
        return [1]


def intraxtbscan(reactants0, entry, reactionpath, count=0):
    scanwd = storage.setting.searching_setting.general.scanwd
    output_d = {}
    uhf_l = check_scan_uhf_intraA(entry)
    struc = entry["struc"]
    s = storage.allstruc.get(struc)
    coordinate = s.struc_info.general.coordinate
    ele_l = s.struc_info.general.ele_l
    at_l = entry["at_combi"]
    push_l = build_push_l(coordinate, ele_l, at_l)
    charge = s.struc_info.general.charge
    nel = s.struc_info.general.nel
    m = s.struc_info.general.multiplicity

    # preparing scan inp
    constrain_setting = {
        "constrain": {
            "force_constant": storage.setting.xtb_setting.k
        },
        "scan": {
            "mode": "concerted" if len(at_l) == 2 else None,
            "distance": push_l if push_l else None,
            "step": storage.setting.xtb_setting.step
        }
    }
    if storage.setting.searching_setting.intra.loose_bond and len(ele_l) > 2:
        # searching weakest bond
        allbonds = SDP([struc],
                       exclude_elements=storage.setting.searching_setting.general.exclude_elements,
                       exclude_pairs=storage.setting.searching_setting.general.exclude_pairs
                       )
        allbonds = sorted(allbonds, key=lambda x: x["mbo"])
        item = allbonds[0]
        if item["mbo"] <= 0.9:
            at1, at2 = item["at_combi"]
            at1, at2 = int(at1), int(at2)
            dd = np.round(Coord(coordinate, ele_l).get_ats_distance(int(at1), int(at2)), 3)
            constrain_setting["constrain"]["loose_bond"] = at1 + 1, at2 + 1, dd

    unreacted = reactants0.copy()
    unreacted.remove(struc)
    unreacted = sorted(unreacted, key=int, reverse=True)
    for u in uhf_l:
        print(f"Simulating intraA path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...",end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")
        create_folder(scan_wd)
        write_xyz(coordinate, ele_l, "coord")
        with open(join(scan_wd, "atl"), "w") as f:
            f.write(json.dumps(entry))

        # call xtb scan
        status = CallxTB(scan_wd,
                         custom_settings={"charge": charge, "uhf": u,
                                          "etemp": storage.setting.xtb_setting.scan_etemp}).scan(
            constrain_setting)
        if status == 0:
            endxyz, ele_l_end = read_coord(join(scan_wd, "xtbopt.xyz"))
            bo = read_bo(join(scan_wd, "wbo"))
            constrain_l = get_topo(endxyz, ele_l_end, bo)
            atcharge = read_atcharge(join(scan_wd, "charges"))

            if storage.setting.searching_setting.general.sep_mode != 0:
                # do a unconstarined opt
                optwd = join(scan_wd, "opt")
                create_folder(optwd)
                shutil.copyfile(join(scan_wd, 'xtbopt.xyz'), join(optwd, 'coord.xyz'))
                status2 = CallxTB(optwd, custom_settings={"charge": charge, "uhf": u}).opt()
                if status2 == 0:
                    endxyz, ele_l_end = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l = get_topo(endxyz, ele_l_end, bo)
                    atcharge = read_atcharge(join(optwd, "charges"))

            output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l,
                                                    "simuhf": u},  # reactants info
                                                   {"atomic_charge": atcharge,
                                                    "charge": charge,
                                                    "nel": nel,
                                                    "uhf": u,
                                                    "coordinate": endxyz,
                                                    "ele_l": ele_l_end,
                                                    "bo": bo,
                                                    "atomic_connection": constrain_l},  # endxyz info
                                                   {"type": "intraA"}  # reaction info
                                                   ]
            print(f"done.")
        else:
            print(f"failed.")
        count += 1
    return output_d, count


def intrarmbff(reactants0, entry, reactionpath, count=0):
    scanwd = storage.setting.searching_setting.general.scanwd
    output_d = {}
    uhf_l = check_scan_uhf_intraA(entry)
    struc = entry["struc"]
    s = storage.allstruc.get(struc)
    coordinate = s.struc_info.general.coordinate
    ele_l = s.struc_info.general.ele_l
    at_l = entry["at_combi"]
    push_l = build_push_l(coordinate, ele_l, at_l)
    charge = s.struc_info.general.charge
    nel = s.struc_info.general.nel
    m = s.struc_info.general.multiplicity

    # get atomic connection (COPY!)
    atom_connection = copy.deepcopy(s.struc_info.general.atomic_connection)
    for item in push_l:
        at1, at2 = int(item[0]), int(item[1])
        value = np.round(float(item[3]), 3)
        atom_connection.append({"atoms": [at1, at2], "value": value})

    unreacted = reactants0.copy()
    unreacted.remove(struc)
    unreacted = sorted(unreacted, key=int, reverse=True)
    for u in uhf_l:
        print(f"Simulating intraA path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...",end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")
        create_folder(scan_wd)

        write_xyz(coordinate, ele_l, "coord")
        with open(join(scan_wd, "atl"), "w") as f:
            f.write(json.dumps(entry))

        with open(join(scan_wd, "newbond"), "w") as f:
            for item in atom_connection:
                f.write(json.dumps(item) + "\n")

        # try first forcing set_angles=Ture to generate target bond connection
        mb_test, mb_ele_test = get_endxyz_from_mb(atom_connection, coordinate, ele_l,
                                                  set_angles=True,
                                                  repulsion_charge=100)
        connectd_trigger = True
        for item in push_l:
            at1, at2 = int(item[0]), int(item[1])
            value = np.round(float(item[3]), 2)
            bond_length = Coord(mb_test, mb_ele_test).get_ats_distance(at1 - 1, at2 - 1)
            if bond_length > value * 1.05:
                connectd_trigger = False
                break
        if not connectd_trigger:
            # otherwise do 2. FF constrain opt
            mbxyz0, ele_l0 = get_endxyz_from_mb(atom_connection, coordinate, ele_l,
                                                set_angles=False,
                                                repulsion_charge=100)
            write_xyz(mbxyz0, ele_l0, "coord_mb0")
            mbxyz0 = pd.read_table(join(scan_wd, "coord_mb0.xyz"), header=None, skiprows=2, sep='\s+')
            mbxyz, mb_ele_l = get_endxyz_from_mb(atom_connection, mbxyz0, ele_l0,
                                                 set_angles=True,
                                                 repulsion_charge=100)
        else:
            mbxyz, mb_ele_l = mb_test, mb_ele_test
        write_xyz(mbxyz, mb_ele_l, "coord_mb")
        status = CallxTB(scan_wd,
                         custom_settings={"charge": charge, "uhf": u, "coordfn": "coord_mb",
                                          "k": storage.setting.xtb_setting.k,
                                          "etemp": storage.setting.xtb_setting.scan_etemp}
                         ).constrained_opt(newbond=atom_connection[-len(push_l):])
        if status == 0:
            endxyz, ele_l_end = read_coord(join(scan_wd, "xtbopt.xyz"))
            bo = read_bo(join(scan_wd, "wbo"))
            constrain_l = get_topo(endxyz, ele_l_end, bo)
            atcharge = read_atcharge(join(scan_wd, "charges"))
            if storage.setting.searching_setting.general.sep_mode != 0:
                # do a unconstarined opt
                optwd = join(scan_wd, "opt")
                create_folder(optwd)
                shutil.copyfile(join(scan_wd, 'xtbopt.xyz'), join(optwd, 'coord.xyz'))
                status2 = CallxTB(optwd, custom_settings={"charge": charge, "uhf": u}).opt()
                if status2 == 0:
                    endxyz, ele_l_end = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l = get_topo(endxyz, ele_l_end, bo)
                    atcharge = read_atcharge(join(optwd, "charges"))

            output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l,
                                                    "simuhf": u},  # reactants info
                                                   {"atomic_charge": atcharge,
                                                    "charge": charge,
                                                    "nel": nel,
                                                    "uhf": u,
                                                    "coordinate": endxyz,
                                                    "ele_l": ele_l_end,
                                                    "bo": bo,
                                                    "atomic_connection": constrain_l},  # endxyz info
                                                   {"type": "intraA"}  # reaction info
                                                   ]
            print(f"done.")
        else:
            print(f"failed.")
        count += 1
    return output_d, count
