# Reaction simulation for intramolecular association
from os.path import join
import shutil
import json
import numpy as np
import pandas as pd
from aRST.script.toolbox import create_folder, write_xyz, read_coord, read_bo, read_atcharge
from aRST.script.spr import SDP
from aRST.script.geom import Coord
from aRST.script.callsys import CallxTB
from aRST.script.do_scan import get_endxyz_from_mb
import aRST.script.initialization as storage
from aRST.script.head import print_simulation


def get_intraD_react_l(reactants0):
    if storage.setting.get_value("searching_setting.disso.at1") is not None:
        # predifined reaction, used in one cycle
        at1 = storage.setting.get_value("searching_setting.disso.at1")
        at2 = storage.setting.get_value("searching_setting.disso.at2")
        react_l = {"struc": [reactants0[0]],
                   "at_combi": [[str(at1), str(at2)]],
                   "gap": 0}
    else:
        react_l = SDP(reactants0,
                      storage.setting.get_value("searching_setting.general.hydrogen", True),
                      storage.setting.get_value("searching_setting.general.halogen", True),
                      ).build_react_l(storage.setting.get_value("searching_setting.disso.mode", 0))
    cut = 2 if storage.setting.get_value("searching_setting.disso.concerted_pair", False) else 1
    for item in react_l:
        pairs = item["at_combi"][:cut]
        item["at_combi"] = pairs
    return react_l


def check_scan_uhf_D(entry):
    struc = entry["struc"]
    nel = storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
    if nel % 2 == 0:
        return [0, 2]
    elif nel % 2 != 0:
        return [1]


def determine_pulling_distance(entry):
    d = storage.setting.get_value("searching_setting.disso.pulling_distance", "auto")
    struc = entry["struc"]
    coord = storage.allstruc.get_value(f"{struc}.struc_info.general.coord")
    at_l = entry["at_combi"]

    push_l = []
    if d == "auto":
        if len(at_l) > 1:
            for pair in at_l[:2]:
                # storage.setting.get_value("searching_setting.disso.pulling_distance", 2) default as 2x covalent radii
                at1, at2 = pair
                d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 2)
                dend = 5 * d12
                # push_l requires atom index starting with 1
                at1 = int(at1) + 1
                at2 = int(at2) + 1
                push_l.append([at1, at2, d12, dend])

        else:
            from aRST.script.toolbox import GraphAnalyzer
            bo = storage.allstruc.get_value(f"{struc}.struc_info.reactivity.bo")
            at1, at2 = at_l[0]
            d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 2)
            # exam topo structure of reactant
            dend = GraphAnalyzer(coord, bo).get_t2t_path(int(at1), int(at2))
            at1 = int(at1) + 1
            at2 = int(at2) + 1
            push_l.append([at1, at2, d12, dend])
    else:
        for pair in at_l[:2]:
            at1, at2 = pair
            d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 2)
            dend = d * d12
            at1 = int(at1) + 1
            at2 = int(at2) + 1
            push_l.append([at1, at2, d12, dend])
    return push_l


def dissoxtbscan(reactants0, entry, reactionpath, count=0):
    output_d = {}
    uhf_l = check_scan_uhf_D(entry)
    for u in uhf_l:
        scan_wd = join(storage.setting.scanwd, f"{reactionpath}_{count}")
        create_folder(scan_wd)

        struc = entry["struc"]
        coord = storage.allstruc.get_value(f"{struc}.struc_info.general.coord")
        write_xyz(coord, "coord")
        at_l = entry["at_combi"]

        # prepare scan
        charge = storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
        nel = storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
        m = storage.allstruc.get_value(f"{struc}.struc_info.general.multiplicity")

        # preparing scan inp
        push_l = determine_pulling_distance(entry)
        constrain_setting = {
            "constrain": {
                "force_constant": storage.setting.get_value("xtb_setting.k", 1)
            },
            "scan": {
                "mode": "concerted" if len(at_l) == 2 else None,
                "distance": push_l if push_l else None,
                "step": storage.setting.get_value("xtb_setting.step", 20)
            }
        }

        status = CallxTB(scan_wd,
                         custom_settings={"charge": charge, "uhf": u,
                                          "etemp": storage.setting.get_value("xtb_setting.scan_etemp", 4000)}).scan(
            constrain_setting)
        if status == 0:
            endxyz = read_coord(join(scan_wd, "xtbopt.xyz"))
            bo = read_bo(join(scan_wd, "wbo"))
            constrain_l, ele_l = CallxTB(scan_wd).get_topo()
            atcharge = read_atcharge(join(scan_wd, "charges"))
            if storage.setting.get_value("searching_setting.general.sep_mode", 0) != 0:
                # do a unconstarined opt
                optwd = join(scan_wd, "opt")
                create_folder(optwd)
                shutil.copyfile(join(scan_wd, 'xtbopt.xyz'), join(optwd, 'coord.xyz'))
                status2 = CallxTB(optwd,
                                  custom_settings={"charge": charge, "uhf": u}).opt()
                if status2 == 0:
                    endxyz = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l, ele_l = CallxTB(optwd).get_topo()
                    atcharge = read_atcharge(join(optwd, "charges"))

            unreacted = reactants0.copy()
            unreacted.remove(struc)
            unreacted = sorted(unreacted, key=int, reverse=True)
            output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l, "simuhf": u},  # reactants info
                                                   {"atomic_charge": atcharge,
                                                    "charge": charge,
                                                    "nel": nel,
                                                    "uhf": u,
                                                    "coord": endxyz,
                                                    "bo": bo,
                                                    "atomic_connection": (constrain_l, ele_l)}  # endxyz info
                                                   ]
            print_simulation(entry["struc"],
                             at_l,
                             [storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
                              for struc in entry["struc"]],
                             [storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
                              for struc in entry["struc"]],
                             reactionpath, count, nel, charge, u
                             )
        count += 1
    return output_d, count


def dissombff(reactants0, entry, reactionpath, count=0):
    output_d = {}
    uhf_l = check_scan_uhf_D(entry)
    for u in uhf_l:
        scan_wd = join(storage.setting.scanwd, f"{reactionpath}_{count}")
        create_folder(scan_wd)

        struc = entry["struc"][0]
        coord = storage.allstruc.get_value(f"{struc}.struc_info.general.coord")
        write_xyz(coord, "coord")
        at_l = entry["at_combi"]

        # prepare scan
        charge = storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
        nel = storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
        m = storage.allstruc.get_value(f"{struc}.struc_info.general.multiplicity")

        # preparing scan inp
        push_l = determine_pulling_distance(entry)
        with open(join(scan_wd, "atl"), "w") as f:
            f.write(json.dumps(entry))
        # get atomic connection
        atom_connection, ele = storage.allstruc.get_value(f"{struc}.struc_info.general.atomic_connection")
        with open(join(scan_wd, "newbond"), "w") as f:
            f.write("")

        for item in push_l:
            at1 = int(item[0])
            at2 = int(item[1])
            value = np.round(float(item[3]), 2)
            atom_connection.append({"atoms": [at1, at2], "value": value})
            with open(join(scan_wd, "newbond"), "a+") as f:
                f.write(json.dumps({"atoms": [at1, at2], "value": value}))
        coordinates = [list(coord.loc[line, 1:4]) for line in range(coord.shape[0])]

        mbxyz0 = get_endxyz_from_mb(atom_connection, coordinates, ele, set_angles=False)
        write_xyz(mbxyz0, "coord_mb0")
        mbxyz0 = pd.read_table(join(scan_wd, "coord_mb0.xyz"), header=None, skiprows=2, sep='\s+')
        coordinates = [list(mbxyz0.loc[line, 1:4]) for line in range(mbxyz0.shape[0])]
        mbxyz = get_endxyz_from_mb(atom_connection, coordinates, ele, set_angles=True)
        write_xyz(mbxyz, "coord_mb")

        status = CallxTB(scan_wd,
                         custom_settings={"charge": charge, "uhf": u, "coordfn": "coord_mb",
                                          "etemp": storage.setting.get_value("xtb_setting.scan_etemp", 4000)
                                          }).constrained_opt(newbond=atom_connection[-len(push_l):])
        if status == 0:
            endxyz = read_coord(join(scan_wd, "xtbopt.xyz"))
            bo = read_bo(join(scan_wd, "wbo"))
            constrain_l, ele_l = CallxTB(scan_wd).get_topo()
            atcharge = read_atcharge(join(scan_wd, "charges"))
            if storage.setting.get_value("searching_setting.general.sep_mode", 0) != 0:
                # do a unconstarined opt
                optwd = join(scan_wd, "opt")
                create_folder(optwd)
                shutil.copyfile(join(scan_wd, 'xtbopt.xyz'), join(optwd, 'coord.xyz'))
                status2 = CallxTB(optwd,
                                  custom_settings={"charge": charge, "uhf": u,
                                                   "etemp": storage.setting.get_value("xtb_setting.scan_etemp",
                                                                                      4000)}).opt()
                if status2 == 0:
                    endxyz = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l, ele_l = CallxTB(optwd).get_topo()
                    atcharge = read_atcharge(join(optwd, "charges"))

            unreacted = reactants0.copy()
            unreacted.remove(struc)
            unreacted = sorted(unreacted, key=int, reverse=True)
            output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l, "simuhf": u},  # reactants info
                                                   {"atomic_charge": atcharge,
                                                    "charge": charge,
                                                    "nel": nel,
                                                    "uhf": u,
                                                    "coord": endxyz,
                                                    "bo": bo,
                                                    "atomic_connection": (constrain_l, ele_l)}  # endxyz info
                                                   ]
            print_simulation(entry["struc"],
                             at_l,
                             [storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
                              for struc in entry["struc"]],
                             [storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
                              for struc in entry["struc"]],
                             reactionpath, count, nel, charge, u
                             )
        count += 1
    return output_d, count
