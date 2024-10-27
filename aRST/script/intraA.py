# Reaction simulation for intramolecular association
import os

import numpy as np
import pandas as pd
from aRST.script.toolbox import create_folder, write_xyz, read_topomol,round_charge
from aRST.script.spr import SPR_HL, build_bo_d

from aRST.script.geom import Coord
from aRST.script.dictionaryInspector import DictionaryInspector


def simulation_interA(reactants,  # reactants list
                      struc_rec,  # correspongding struc info in dict
                      searching_setting,
                      xtb_setting,
                      wd,  # working folder
                      reactionpath,  # str numbering for prior path
                      count=0):  # :return: output_d[count] = [{reactants info},{endxyz info}]

    output_d = {}

    struc_inspector = DictionaryInspector(struc_rec)
    ###### rangking reactivity ######
    react_aa, react_bb, react_ab, react_ba = SPR_HL(reactants,
                                                    struc_rec,
                                                    allow_hydrogen=searching_setting.get("hydrogen", True),
                                                    allow_halogen=searching_setting.get("halogen", True),
                                                    search_range=0
                                                    )  # react_d:  key= 5.93, value = {('1', '1'), [['0', '0']]}

    bo = build_bo_d(reactants,
                    struc_inspector,
                    allow_hydrogen=searching_setting.get("hydrogen"),
                    allow_halogen=searching_setting.get("halogen", True))  # bond connection info
    if not bo:
        return output_d

    intraa = {}
    for react_d in [react_aa, react_bb, react_ab, react_ba]:
        for current_key, current_item in react_d.items():
            #  5.93 {('1', '1'): [['0', '0']]}
            for struc_combi, at_l in current_item.items():
                struc = struc_combi[0]
                if 1 not in [len(set(k)) for k in at_l]:  # not same atom
                    for pair in at_l:
                        if pair in bo[struc] or pair[::-1] in bo[struc]:  # not bonded
                            continue
                        else:
                            intraa[current_key] = struc, at_l

    # react_d:  key= 5.93, value = {('1', '1'), [['0', '0']]}
    # intraa:   key= 5.93,value = ('1',[['0', '0']])

    if not intraa:
        raise Exception(f"no possible intraA possibilities can be found in reactants list: {reactants}.")

    scanmode = searching_setting.get("scan_mode", 0)
    nintra = searching_setting.get("number", 1)
    from aRST.script.callsys import CallxTB
    from aRST.script.do_scan import build_push_l
    if scanmode == 0:
        # standard xtbscan
        for gap in sorted(intraa.keys())[:nintra]:
            scan_wd = os.path.join(wd, "scan", f"{reactionpath}_{count}")
            create_folder(scan_wd)

            strucid = intraa[gap][0]
            coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            write_xyz(coord, f"coord")

            # prepare scan setteing
            charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            nel = struc_inspector.get_value(f"{strucid}.struc_info.general.nel")
            uhf = struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1

            at_l = intraa[gap][1]
            push_l = build_push_l(coord, at_l)

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
                wboe = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.bo')
                tmpl = [i if i > 0.3 else np.inf for i in wboe.loc[:, 2].tolist()]
                minbo = min(tmpl)
                if 0.3 < minbo < 0.9:
                    at1, at2 = wboe.loc[tmpl.index(minbo), 0], wboe.loc[tmpl.index(minbo), 1]
                    dd = np.round(Coord(coord).get_ats_distance(int(at1 - 1), int(at2 - 1)), 3)
                    constrain_setting["constrain"]["loose_bond"] = at1, at2, dd

            # call xtb scan
            status = CallxTB(scan_wd, {"method": xtb_setting.get("method", "gfn2"),
                                       "coordfn": "coord",
                                       "optlevel": xtb_setting.get("optlevel", None),
                                       "charge": charge,
                                       "uhf": uhf,
                                       "alpb": xtb_setting.get("alpb", "h2o"),
                                       "etemp": xtb_setting.get("scan_etemp", "4000")
                                       }).scan(constrain_setting)

            if status == 0:
                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                constrain_l, ele = read_topomol(os.path.join(scan_wd, 'xtbtopo.mol')) if os.path.exists(
                            os.path.join(scan_wd, 'xtbtopo.mol')) else [],[]
                if os.path.getsize(os.path.join(scan_wd, "wbo")) != 0:
                    bo = pd.read_csv(os.path.join(scan_wd, "wbo"), header=None, sep='\s+')
                else:
                    bo = None
                atcharge = pd.read_table(os.path.join(scan_wd, "charges"), header=None, sep='\s+')
                output_d[f"{count}"] = [{"reactants": [strucid],
                                         "unreacted": [st for st in reactants if st != strucid]
                                         if len(reactants) > 1 else []},  # reactants info
                                        {"atomic_charge": atcharge,
                                         "charge":round_charge(sum(atcharge),0.5),
                                         "nel": nel,
                                         "coord": endxyz,
                                         "bo": bo,
                                         "atomic_connection": (constrain_l, ele)}  # endxyz info
                                        ]
            count += 1
        return output_d
    elif scanmode == 1:  # use molbar FF
        from aRST.script.do_scan import get_endxyz_from_mb

        for gap in sorted(intraa.keys())[:nintra]:
            scan_wd = os.path.join(wd, "scan", f"{reactionpath}_{count}")
            create_folder(scan_wd)

            strucid = intraa[gap][0]
            coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            write_xyz(coord, f"coord")

            # prepare scan setteing
            charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            nel = struc_inspector.get_value(f"{strucid}.struc_info.general.nel")
            uhf = struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1

            at_l = intraa[gap][1]
            push_l = build_push_l(coord, at_l)

            # get atomic connection
            atom_connection, ele = struc_inspector.get_value(f"{strucid}.struc_info.atomic_connection")
            for item in push_l:
                at1, at2 = int(item[0]), int(item[1])
                value = float(item[3])
                atom_connection.append({"atoms": [at1, at2], "value": value})
            coordinates = [list(coord.loc[line, 1:4]) for line in range(coord.shape[0])]
            mbxyz = get_endxyz_from_mb(atom_connection, coordinates, ele)
            write_xyz(mbxyz, "coord0")
            status = CallxTB(scan_wd, {"method": xtb_setting.get("method", "gfn2"),
                                       "coordfn": "coord0",
                                       "optlevel": xtb_setting.get("optlevel", None),
                                       "charge": charge,
                                       "uhf": uhf,
                                       "alpb": xtb_setting.get("alpb", "h2o"),
                                       "etemp": xtb_setting.get("scan_etemp", "4000")
                                       }).opt()
            if status == 0:
                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                constrain_l, ele = read_topomol(os.path.join(scan_wd, 'xtbtopo.mol')) if os.path.exists(
                            os.path.join(scan_wd, 'xtbtopo.mol')) else [],[]
                if os.path.getsize(os.path.join(scan_wd, "wbo")) != 0:
                    bo = pd.read_csv(os.path.join(scan_wd, "wbo"), header=None, sep='\s+')
                else:
                    bo = None
                atcharge = pd.read_table(os.path.join(scan_wd, "charges"), header=None, sep='\s+')
                output_d[f"{count}"] = [{"reactants": [strucid],
                                         "unreacted": [st for st in reactants if st != strucid]
                                         if len(reactants) > 1 else []},  # reactants info
                                        {"atomic_charge": atcharge,
                                         "charge":round_charge(sum(atcharge),0.5),
                                         "nel": nel,
                                         "coord": endxyz,
                                         "bo": bo,
                                         "atomic_connection": (constrain_l, ele)}  # endxyz info
                                        ]
            count += 1
        return output_d
    else:
        raise Exception("unknown keyword for scan_mode, accecpted values include (0,1).")
