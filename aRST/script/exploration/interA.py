# Reaction simulation for intramolecular association
import json
import sys
from os.path import join
import shutil
import copy
import numpy as np
import pandas as pd
from aRST.script.toolbox import create_folder, write_xyz, read_coord, read_bo, read_atcharge
from aRST.script.spr import SPR, SDP
from aRST.script.geom import Coord, alignment
from aRST.script.callsys import CallxTB
from aRST.script.do_scan import build_push_l
from aRST.script.do_scan import get_endxyz_from_mb
import aRST.script.initialization as storage
from aRST.script.head import print_simulation, print_react_l


def get_interA_react_l(reactants0,printall=True):
    react_l = SPR(reactants0,
                  storage.setting.get_value("searching_setting.general.hydrogen", True),
                  storage.setting.get_value("searching_setting.general.halogen", True)
                  ).inter(storage.setting.get_value("searching_setting.general.spr_mode", 0))
    cut = 2 if storage.setting.get_value("searching_setting.inter.concerted_pair", False) else 1
    for item in react_l:
        pairs = item["at_combi"][:cut]
        item["at_combi"] = pairs

    # print top 20 react_l
    if printall:
        print_react_l(reactants0, react_l, rtype="interA")
        # write full react_l in file
        if storage.setting.get_value("searching_setting.general.general", 1) == 1:
            with open(join(storage.setting.wd0, "react_l"), "w") as f:
                for item in react_l:
                    f.write(json.dumps(item) + "\n")

    if storage.setting.get_value("searching_setting.inter.at1") is not None:
        # overwrite react_l with user's input
        # predifined reaction, used in one cycle
        at1 = storage.setting.get_value("searching_setting.inter.at1")
        at2 = storage.setting.get_value("searching_setting.inter.at2")
        if storage.setting.get_value("searching_setting.inter.at3") is not None:
            at3 = storage.setting.get_value("searching_setting.inter.at3")
            at4 = storage.setting.get_value("searching_setting.inter.at4")
            react_l = [{"struc_combi": [reactants0[0], reactants0[1]],
                       "at_combi": [[str(at1), str(at2)],[str(at3), str(at4)]],
                       "gap": 0}]
            storage.setting.set_value("searching_setting.inter.concerted_pair", True)
        else:
            react_l = [{"struc_combi": [reactants0[0], reactants0[1]],
                   "at_combi": [[str(at1), str(at2)]],
                   "gap": 0}]

    return react_l


def generate_reactive_complex(entry, reactionpath, count=0):
    # prepare scan_wd
    scan_wd = join(storage.setting.scanwd, f"{reactionpath}_{count}")
    create_folder(scan_wd)

    # prepare scan_xyz

    struc1, struc2 = entry["struc_combi"]
    coord1 = storage.allstruc.get_value(f"{struc1}.struc_info.general.coord")
    coord2 = storage.allstruc.get_value(f"{struc2}.struc_info.general.coord")
    nat_coord1 = int(coord1.shape[0])

    # generate reactiv complex
    at_l = entry["at_combi"]
    at1, at2 = at_l[0]
    at1 = int(at1)
    at2 = int(at2)
    alignedcoord = alignment(coord1, coord2, at1, at2,
                             extra_setting={"maxovlap": storage.setting.get_value(
                                 "searching_setting.inter.maxovlap", False),
                                 "aligning_distance": storage.setting.get_value(
                                     "searching_setting.inter.aligning_distance", 5)})
    write_xyz(alignedcoord, f"coord")

    # prepare scan setteing
    charge = sum([storage.allstruc.get_value(f"{k}.struc_info.general.charge") for k in [struc1, struc2]])
    nelA = storage.allstruc.get_value(f"{struc1}.struc_info.general.nel")
    nelB = storage.allstruc.get_value(f"{struc2}.struc_info.general.nel")

    push_l = build_push_l(alignedcoord, at_l, nat_coord1)

    return scan_wd, struc1, struc2, at_l, alignedcoord, nat_coord1, charge, nelA + nelB, push_l


def check_scan_uhf_interA(entry):
    struc1, struc2 = entry["struc_combi"]
    nelA = storage.allstruc.get_value(f"{struc1}.struc_info.general.nel")
    nelB = storage.allstruc.get_value(f"{struc2}.struc_info.general.nel")
    if nelA % 2 == 0 and nelB % 2 == 0:
        return [0]
    elif (nelA + nelB) % 2 != 0:
        return [1]
    else:
        return [0, 2]


def interxtbscan(reactants0, entry, reactionpath, count=0):
    output_d = {}
    # preset output
    uhf_l = check_scan_uhf_interA(entry)
    for u in uhf_l:
        scan_wd, \
        struc1, struc2, at_l, \
        alignedcoord, nat_coord1, charge, totnel, push_l = generate_reactive_complex(entry, reactionpath, count=count)

        with open(join(scan_wd, "atl"), "w") as f:
            f.write(json.dumps(entry))
        # preparing scan inp
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
        if storage.setting.get_value("searching_setting.inter.loose_bond", False):
            # searching weakest bond
            allbonds = SDP([struc1, struc2],
                           storage.setting.get_value("searching_setting.general.hydrogen", True),
                           storage.setting.get_value("searching_setting.general.halogen", True),
                           ).allbonds
            allbonds = sorted(allbonds, key=lambda x: x["mbo"])
            item = allbonds[0]
            if item["mbo"] <= 0.9:
                if item["struc"] == struc1:
                    at1, at2 = item["at_combi"]
                    at1, at2 = int(at1), int(at2)
                else:
                    at1, at2 = item["at_combi"]
                    at1, at2 = int(at1) + nat_coord1, int(at2) + nat_coord1
                dd = np.round(Coord(alignedcoord).get_ats_distance(int(at1), int(at2)), 3)
                constrain_setting["constrain"]["loose_bond"] = at1 + 1, at2 + 1, dd

        # call xtb scan
        status = CallxTB(scan_wd,
                         custom_settings={"charge": charge, "uhf": u, "etemp": storage.setting.get_value(
                             "xtb_setting.scan_etemp", 4000)}).scan(constrain_setting)

        if status == 0:
            endxyz = read_coord(join(scan_wd, "xtbopt.xyz"))
            bo = read_bo(join(scan_wd, "wbo"))
            constrain_l, ele_l = CallxTB(scan_wd).get_topo()
            atcharge = read_atcharge(join(scan_wd, "charges"))
            if storage.setting.get_value("searching_setting.general.sep_mode", 0) != 0:
                # do a unconstarined opt
                optwd = join(scan_wd, "opt")
                create_folder(optwd)
                write_xyz(endxyz, "coord")
                status2 = CallxTB(optwd, custom_settings={"charge": charge, "uhf": u}).opt()
                if status2 == 0:
                    endxyz = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l, ele_l = CallxTB(optwd).get_topo()
                    atcharge = read_atcharge(join(optwd, "charges"))
            # md optional
            if storage.setting.get_value("searching_setting.general.MD", False):
                mdwd = join(scan_wd, "md")
                bo, constrain_l, ele_l, atcharge, endxyz = _short_md_for_product_struc_guess(mdwd, endxyz, push_l,
                                                                                             charge, u)

            unreacted = reactants0.copy()
            unreacted.remove(struc1)
            unreacted.remove(struc2)
            unreacted = sorted(unreacted, key=int, reverse=True)
            output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc1, struc2],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l, "simuhf": u},  # reactants info
                                                   {"atomic_charge": atcharge,
                                                    "charge": charge,
                                                    "nel": totnel,
                                                    "uhf": u,
                                                    "coord": endxyz,
                                                    "bo": bo,
                                                    "atomic_connection": (constrain_l, ele_l)}  # endxyz info
                                                   ]
            print_simulation(entry["struc_combi"],
                             at_l,
                             [storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
                              for struc in entry["struc_combi"]],
                             [storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
                              for struc in entry["struc_combi"]],
                             reactionpath, count, totnel, charge, u, nat_1=nat_coord1
                             )
        count += 1
    return output_d, count


def intermbff(reactants0, entry, reactionpath, count=0):
    output_d = {}
    uhf_l = check_scan_uhf_interA(entry)
    for u in uhf_l:
        scan_wd, \
        struc1, struc2, at_l, \
        alignedcoord, nat_coord1, charge, totnel, push_l = generate_reactive_complex(entry, reactionpath, count=count)

        with open(join(scan_wd, "atl"), "w") as f:
            f.write(json.dumps(entry))
        # get atomic connection
        atom_connection1, ele1 = storage.allstruc.get_value(f"{struc1}.struc_info.general.atomic_connection")
        atom_connection2, ele2 = storage.allstruc.get_value(f"{struc2}.struc_info.general.atomic_connection")
        add_atom_connection2 = copy.deepcopy(atom_connection2)
        for item in add_atom_connection2:
            for at in range(len(item["atoms"])):
                item["atoms"][at] += nat_coord1

        atom_connection = atom_connection1 + add_atom_connection2
        ele = ele1 + ele2
        with open(join(scan_wd, "newbond"), "w") as f:
            f.write("")

        for item in push_l:
            at1, at2 = int(item[0]), int(item[1])
            value = np.round(float(item[3]), 2)
            atom_connection.append({"atoms": [at1, at2], "value": value})
            with open(join(scan_wd, "newbond"), "a+") as f:
                f.write(json.dumps({"atoms": [at1, at2], "value": value}))
        coordinates = [list(alignedcoord.loc[line, 1:4]) for line in range(alignedcoord.shape[0])]
        # mbxyz=get_endxyz_from_mb(atom_connection, coordinates, ele, set_angles=True)
        mbxyz0 = get_endxyz_from_mb(atom_connection, coordinates, ele, set_angles=False)
        write_xyz(mbxyz0, "coord_mb0")
        mbxyz0 = pd.read_table(join(scan_wd, "coord_mb0.xyz"), header=None, skiprows=2, sep='\s+')
        coordinates = [list(mbxyz0.loc[line, 1:4]) for line in range(mbxyz0.shape[0])]
        mbxyz = get_endxyz_from_mb(atom_connection, coordinates, ele, set_angles=True)
        write_xyz(mbxyz, "coord_mb")
        # constrained opt in gfn2
        status = CallxTB(scan_wd,
                         custom_settings={"charge": charge, "uhf": u,
                                          "k": storage.setting.get_value("xtb_setting.k", 1),
                                          "coordfn": "coord_mb",
                                          "etemp": storage.setting.get_value(
                                              "xtb_setting.scan_etemp", 4000)
                                          }).constrained_opt(
            newbond=atom_connection[-len(push_l):])
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
                status2 = CallxTB(optwd, custom_settings={"charge": charge, "uhf": u}).opt()
                if status2 == 0:
                    endxyz = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l, ele_l = CallxTB(optwd).get_topo()
                    atcharge = read_atcharge(join(optwd, "charges"))
            if storage.setting.get_value("searching_setting.general.MD", False):
                mdwd = join(scan_wd, "md")
                bo, constrain_l, ele_l, atcharge, endxyz = _short_md_for_product_struc_guess(mdwd, endxyz, push_l,
                                                                                             charge, u)

            unreacted = reactants0.copy()
            unreacted.remove(struc1)
            unreacted.remove(struc2)
            output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc1, struc2],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l,
                                                    "simuhf": u},  # reactants info
                                                   {"atomic_charge": atcharge,
                                                    "charge": charge,
                                                    "nel": totnel,
                                                    "uhf": u,
                                                    "coord": endxyz,
                                                    "bo": bo,
                                                    "atomic_connection": (constrain_l, ele)}  # endxyz info
                                                   ]
            print_simulation(entry["struc_combi"],
                             at_l,
                             [storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
                              for struc in entry["struc_combi"]],
                             [storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
                              for struc in entry["struc_combi"]],
                             reactionpath, count, totnel, charge, u, nat_1=nat_coord1
                             )
        count += 1
    return output_d, count


def _short_md_for_product_struc_guess(mdwd, coord, constrained_bonds, charge, uhf):
    create_folder(mdwd)
    write_xyz(coord, "coord")
    if storage.setting.get_value("searching_setting.md.temp") is not None:
        mdtemp = storage.setting.get_value("searching_setting.md.temp")
    else:
        mdtemp = storage.setting.get_value("xtb_setting.scan_etemp")
    xtbmdparse = CallxTB(mdwd, custom_settings={"charge": charge, "uhf": uhf,
                                                "mdetemp": mdtemp,
                                                "mdtime":storage.setting.get_value("searching_setting.md.time",10),
                                                "mdstep": storage.setting.get_value("searching_setting.md.step", 1),
                                                "mdshake": storage.setting.get_value("searching_setting.md.shake", 1),
                                                "mddump": storage.setting.get_value("searching_setting.md.dump", 50),
                                                "etemp":storage.setting.get_value("xtb_setting.scan_etemp") # sp after md
                                                })
    status = xtbmdparse.short_md(constrained_bonds=constrained_bonds)
    if status != 0:
        raise ValueError(f"md has error, check {mdwd}")
    # get energy min struc in the md simulation
    # coord can't be single at in this case
    coord, emd = xtbmdparse.get_Emin_geom_from_mdtrj(nat=coord.shape[0])
    write_xyz(coord, "emincoord")
    # doing a sp in the same folder for sp information
    status = xtbmdparse.sp()
    if status != 0:
        raise ValueError(f"sp after md has error, check {mdwd}")
    bo = read_bo(join(mdwd, "wbo"))
    constrain_l, ele_l = xtbmdparse.get_topo()
    atcharge = read_atcharge(join(mdwd, "charges"))
    return bo, constrain_l, ele_l, atcharge, coord
