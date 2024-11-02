# Reaction simulation for intramolecular association
import os
import shutil

import numpy as np
import pandas as pd
from aRST.script.toolbox import create_folder, write_xyz, read_topomol, round_charge
from aRST.script.spr import SPR_HL
from aRST.script.dictionaryInspector import DictionaryInspector
from aRST.script.geom import Coord, alignment


def simulation_interA(reactants,  # reactants list
                      struc_rec,  # correspongding struc info in dict
                      searching_setting,
                      xtb_setting,
                      wd,  # working folder
                      reactionpath,  # str numbering for prior path
                      count=0):  # :return: output_d[count] = [{reactants info},{endxyz info}]

    struc_inspector = DictionaryInspector(struc_rec)
    ###### rangking reactivity ######
    react_aa, react_bb, react_ab, react_ba = SPR_HL(reactants,
                                                    struc_rec,
                                                    allow_hydrogen=searching_setting.get("hydrogen", True),
                                                    allow_halogen=searching_setting.get("halogen", True),
                                                    search_mode=searching_setting.get("search_mode", 1),
                                                    )
    # react_d:  key= 5.93, value = {('1', '1'), [['0', '0']]}

    # merge a reaction dict
    intera = {}  # intera:   key= 5.93, value = ('1','1',[['0', '0']])
    for react_d in [react_aa, react_bb, react_ab, react_ba]:
        for current_key, current_item in react_d.items():
            for struc_combi, at_l in current_item.items():
                struc1, struc2 = struc_combi
                intera[current_key] = struc1, struc2, at_l
    if not intera:
        raise Exception(f"no possible interA possibilities can be found in reactants list: {reactants}.")

    # resort intera to screen same scan
    def screen_identical_scan(input_d):
        from collections import Counter
        # before: struc1: provide hoao struc2: provide luao
        # after: struc1 index <struc2 index
        out_d = {}
        for gap in sorted(input_d.keys()):
            if gap not in out_d.keys():
                tmps1, tmps2 = input_d[gap][:2]
                if int(tmps1) <= int(tmps2):
                    s1 = tmps1
                    s2 = tmps2
                    atl = input_d[gap][2]
                else:
                    s1 = tmps2
                    s2 = tmps1
                    atl = []
                    for item in input_d[gap][2]:
                        if not atl:
                            atl.append(list(reversed(item)))
                        elif int(list(reversed(item))[0]) <=int(atl[0][0]):
                            atl.insert(0,list(reversed(item)))
                        else:atl.append(list(reversed(item)))

                trigger = True
                for item in out_d.values():
                    if len(item[2]) !=len(atl):
                        continue
                    else:
                        if item[2][0] == atl[0]:
                            trigger=False
                            break

                if trigger:
                    out_d[gap] = s1, s2, atl
        return out_d

    intera = screen_identical_scan(intera)

    output_d = {}

    ###### perform simulation ######
    # prepare for reaction simulation
    def generate_reactive_complex(gap):
        # prepare scan_wd
        scan_wd = os.path.join(wd, "scan", f"{reactionpath}_{count}")
        create_folder(scan_wd)

        # prepare scan_xyz
        struc1, struc2 = intera[gap][0], intera[gap][1]
        coord1 = struc_inspector.get_value(f"{struc1}.struc_info.general.coord")
        coord2 = struc_inspector.get_value(f"{struc2}.struc_info.general.coord")
        nat_coord1 = int(coord1.shape[0])

        # generate reactiv complex
        at_l = intera[gap][2]
        at1, at2 = int(at_l[0][0]), int(at_l[0][1])
        alignedcoord = alignment(coord1, coord2, at1, at2,
                                 extra_setting={"maxovlap": searching_setting.get("maxovlap", False),
                                                "aligning_distance": searching_setting.get("aligning_distance", 5)})
        write_xyz(alignedcoord, f"coord")

        # prepare scan setteing
        charge = sum([struc_inspector.get_value(f"{k}.struc_info.general.charge") for k in [struc1, struc2]])
        nelA = struc_inspector.get_value(f"{struc1}.struc_info.general.nel")
        nelB = struc_inspector.get_value(f"{struc2}.struc_info.general.nel")
        if nelA % 2 == 0 and nelB % 2 == 0:
            uhf = 0
        elif nelA % 2 != 0 and nelB % 2 != 0:
            uhf = np.inf  # try both 0 and 2 (additional reaction path)
        else:
            uhf = 1

        push_l = build_push_l(alignedcoord, at_l, nat_coord1)

        return scan_wd, struc1, struc2, at_l, alignedcoord, nat_coord1, charge, uhf, nelA + nelB, push_l

    scanmode = searching_setting.get("scan_mode", 0)
    ninter = searching_setting.get("number", 1)

    #########  tmp hacked in  #########
    if ninter == "all":
        ninter = len(intera.keys())

    #########  tmp hacked out  #########

    from aRST.script.callsys import CallxTB

    if scanmode == 0:
        # standard xtbscan
        from aRST.script.do_scan import build_push_l
        for gap in sorted(intera.keys())[:ninter]:
            scan_wd, \
            struc1, struc2, at_l, \
            alignedcoord, nat_coord1, charge, uhf, totnel, push_l = generate_reactive_complex(gap)

            #########  tmp hacked in  #########
            at_l = [at_l[0]]
            with open(os.path.join(scan_wd, "atl"), "w") as f:
                f.write("at1,at2,hlgap(eV)\n")
                f.write("%d,%d,%.4f\n" % (int(at_l[0][0]), int(at_l[0][1]), gap))
            #########  tmp hacked out  #########

            # preparing scan inp
            constrain_setting = {
                "constrain": {
                    "force_constant": xtb_setting.get("k", 1)
                },
                "scan": {
                    "mode": "concerted" if len(at_l) == 2 else None,
                    "distance": push_l if push_l else None,
                    "step": xtb_setting.get("step", 20)
                }
            }
            if searching_setting.get("loose_bond", False):
                # searching weakest bond
                at1, at2, bo = None, None, None
                for struc in [struc1, struc2]:
                    wboe = struc_inspector.get_value(f'{struc}.struc_info.reactivity.bo')
                    tmpl = [i if i > 0.3 else np.inf for i in wboe.loc[:, 2].tolist()]
                    minbo = min(tmpl)
                    if minbo > 0.9:
                        continue
                    elif bo is None or minbo < bo:
                        if struc == struc1:
                            at1, at2 = wboe.loc[tmpl.index(minbo), 0], wboe.loc[tmpl.index(minbo), 1]
                        else:
                            at1, at2 = wboe.loc[tmpl.index(minbo), 0] + nat_coord1, \
                                       wboe.loc[tmpl.index(minbo), 1] + nat_coord1
                if at1 is not None and at2 is not None:
                    dd = np.round(Coord(alignedcoord).get_ats_distance(int(at1 - 1), int(at2 - 1)), 3)
                    constrain_setting["constrain"]["loose_bond"] = at1, at2, dd

            # call xtb scan
            status = CallxTB(scan_wd, {"method": xtb_setting.get("method", "gfn2"),
                                       "coordfn": "coord",
                                       "optlevel": xtb_setting.get("optlevel", None),
                                       "charge": charge,
                                       "uhf": uhf if not np.isinf(uhf) else 0,
                                       "alpb": xtb_setting.get("alpb", "h2o"),
                                       "etemp": xtb_setting.get("scan_etemp", "4000")
                                       }).scan(constrain_setting)

            if status == 0:
                endwd = scan_wd
                if searching_setting.get("sep_mode", 0) !=0:
                    # do a unconstarined opt
                    optwd = os.path.join(scan_wd, "opt")
                    create_folder(optwd)
                    shutil.copyfile(os.path.join(scan_wd, 'xtbopt.xyz'), os.path.join(optwd, 'coord.xyz'))
                    status2 = CallxTB(optwd, {"method": xtb_setting.get("method", "gfn2"),
                                              "coordfn": "coord",
                                              "optlevel": xtb_setting.get("optlevel", None),
                                              "charge": charge,
                                              "uhf": uhf if not np.isinf(uhf) else 0,
                                              "alpb": xtb_setting.get("alpb", "h2o"),
                                              "etemp": xtb_setting.get("struc_etemp", "4000")
                                              }).opt()
                    if status2 == 0:
                        endwd = optwd
                endxyz = pd.read_table(os.path.join(endwd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                constrain_l, ele = read_topomol(os.path.join(endwd, 'xtbtopo.mol')) if os.path.exists(
                    os.path.join(endwd, 'xtbtopo.mol')) else [], []
                bo = pd.read_csv(os.path.join(endwd, "wbo"), header=None, sep='\s+') if os.path.getsize(
                    os.path.join(endwd, "wbo")) != 0 else None
                atcharge = pd.read_table(os.path.join(endwd, "charges"), header=None, sep='\s+')

                output_d[f"{count}"] = [{"reactants": [struc1, struc2],
                                         "unreacted": [st for st in reactants
                                                       if st != struc1 and st != struc2]
                                         if len(reactants) > 2 else []},  # reactants info
                                        {"atomic_charge": atcharge,
                                         "charge": round_charge(sum(atcharge), 0.5),
                                         "nel": totnel,
                                         "coord": endxyz,
                                         "bo": bo,
                                         "atomic_connection": (constrain_l, ele)}  # endxyz info
                                        ]
            count += 1

            if np.isinf(uhf):
                shutil.copytree(os.path.join(wd, "scan", f"{reactionpath}_{count - 1}"),
                                os.path.join(wd, "scan", f"{reactionpath}_{count}"))
                scan_wd = os.path.join(wd, "scan", f"{reactionpath}_{count}")
                os.chdir(scan_wd)

                status = CallxTB(scan_wd, {"method": xtb_setting.get("method", "gfn2"),
                                           "coordfn": "coord",
                                           "optlevel": xtb_setting.get("optlevel", None),
                                           "charge": charge,
                                           "uhf": 2,
                                           "alpb": xtb_setting.get("alpb", "h2o"),
                                           "etemp": xtb_setting.get("scan_etemp", "4000")
                                           }).scan(constrain_setting)
                if status == 0:
                    endwd = scan_wd
                    if searching_setting.get("sep_mode", 0) != 0:
                        # do a unconstarined opt
                        optwd = os.path.join(scan_wd, "opt")
                        create_folder(optwd)
                        shutil.copyfile(os.path.join(scan_wd, 'xtbopt.xyz'), os.path.join(optwd, 'coord.xyz'))
                        status2 = CallxTB(optwd, {"method": xtb_setting.get("method", "gfn2"),
                                                  "coordfn": "coord",
                                                  "optlevel": xtb_setting.get("optlevel", None),
                                                  "charge": charge,
                                                  "uhf": 2,
                                                  "alpb": xtb_setting.get("alpb", "h2o"),
                                                  "etemp": xtb_setting.get("scan_etemp", "4000")
                                                  }).opt()
                        if status2 == 0:
                            endwd = optwd
                    endxyz = pd.read_table(os.path.join(endwd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                    constrain_l, ele = read_topomol(os.path.join(endwd, 'xtbtopo.mol')) if os.path.exists(
                        os.path.join(endwd, 'xtbtopo.mol')) else [], []
                    bo = pd.read_csv(os.path.join(endwd, "wbo"), header=None, sep='\s+') if os.path.getsize(
                        os.path.join(endwd, "wbo")) != 0 else None
                    atcharge = pd.read_table(os.path.join(endwd, "charges"), header=None, sep='\s+')

                    output_d[f"{count}"] = [{"reactants": [struc1, struc2],
                                             "unreacted": [st for st in reactants
                                                           if st != struc1 and st != struc2]
                                             if len(reactants) > 2 else []},  # reactants info
                                            {"atomic_charge": atcharge,
                                             "charge": round_charge(sum(atcharge), 0.5),
                                             "nel": totnel,
                                             "coord": endxyz,
                                             "bo": bo,
                                             "atomic_connection": (constrain_l, ele)}  # endxyz info
                                            ]
                count += 1
        return output_d
    elif scanmode == 1:  # use molbar FF
        from aRST.script.do_scan import get_endxyz_from_mb

        for gap in sorted(intera.keys())[:ninter]:
            scan_wd, \
            struc1, struc2, at_l, \
            alignedcoord, nat_coord1, charge, uhf, totnel, push_l = generate_reactive_complex(gap)

            # get atomic connection
            atom_connection1, ele1 = struc_inspector.get_value(f"{struc1}.struc_info.atomic_connection")
            atom_connection2, ele2 = struc_inspector.get_value(f"{struc2}.struc_info.atomic_connection")
            for item in atom_connection2:
                for at in item["atoms"]:
                    at += nat_coord1

            atom_connection = atom_connection1 + atom_connection2
            ele = ele1 + ele2
            for item in push_l:
                at1, at2 = int(item[0]), int(item[1])
                value = float(item[3])
                atom_connection.append({"atoms": [at1, at2], "value": value})
            coordinates = [list(alignedcoord.loc[line, 1:4]) for line in range(alignedcoord.shape[0])]
            mbxyz = get_endxyz_from_mb(atom_connection, coordinates, ele)
            write_xyz(mbxyz, "coord0")
            status = CallxTB(scan_wd, {"method": xtb_setting.get("method", "gfn2"),
                                       "coordfn": "coord0",
                                       "optlevel": xtb_setting.get("optlevel", None),
                                       "charge": charge,
                                       "uhf": 1 if totnel % 2 != 0 else 2,
                                       "alpb": xtb_setting.get("alpb", "h2o"),
                                       "etemp": xtb_setting.get("scan_etemp", "4000")
                                       }).opt()
            if status == 0:
                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                constrain_l, ele = read_topomol(os.path.join(scan_wd, 'xtbtopo.mol')) if os.path.exists(
                    os.path.join(scan_wd, 'xtbtopo.mol')) else [], []
                if os.path.getsize(os.path.join(scan_wd, "wbo")) != 0:
                    bo = pd.read_csv(os.path.join(scan_wd, "wbo"), header=None, sep='\s+')
                else:
                    bo = None
                atcharge = pd.read_table(os.path.join(scan_wd, "charges"), header=None, sep='\s+')
                output_d[f"{count}"] = [{"reactants": [struc1, struc2],
                                         "unreacted": [st for st in reactants
                                                       if st != struc1 and st != struc2]
                                         if len(reactants) > 2 else []},  # reactants info
                                        {"atomic_charge": atcharge,
                                         "charge": round_charge(sum(atcharge), 0.5),
                                         "nel": totnel,
                                         "coord": endxyz,
                                         "bo": bo,
                                         "atomic_connection": (constrain_l, ele)}  # endxyz info
                                        ]
            count += 1
        return output_d
    else:
        raise Exception("unknown keyword for scan_mode, accecpted values include (0,1).")
