# Reaction simulation for intramolecular association
import os
import shutil

import numpy as np
import pandas as pd
from aRST.script.toolbox import create_folder, write_xyz, read_topomol,round_charge
from aRST.script.spr import SPR_dissociation
from aRST.script.dictionaryInspector import DictionaryInspector
from aRST.script.geom import Coord


def simulation_intraD(reactants,  # reactants list
                      struc_rec,  # correspongding struc info in dict
                      searching_setting,
                      xtb_setting,
                      wd,  # working folder
                      reactionpath,  # str numbering for prior path
                      count=0):  # :return: output_d[count] = [{reactants info},{endxyz info}]


    struc_inspector = DictionaryInspector(struc_rec)
    output_d = {}

    ###### rangking reactivity ######
    BO = SPR_dissociation(reactants,  # reactants id list
                          struc_rec,  # struc dict record
                          allow_hydrogen=searching_setting.get("hydrogen", True),
                          allow_halogen=searching_setting.get("halogen", True),
                          search_mode=0 if str.lower(searching_setting.get("mode", "wboe")) == "wboe" else 1)
    if not BO:
        raise Exception(f"no possible dissociation possibilities can be found in reactants list: {reactants}.")

    scanmode = searching_setting.get("scan_mode", 0)
    ndisso = searching_setting.get("number", 1)

    if scanmode == 0:
        # standard xtbscan
        from aRST.script.callsys import CallxTB
        for dgap in sorted(BO.keys())[:ndisso]:
            # prepare scan_wd
            scan_wd = os.path.join(wd, "scan", f"{reactionpath}_{count}")
            create_folder(scan_wd)

            # prepare scan_xyz
            info_l = BO[dgap][0]
            strucid = info_l[0]
            coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            write_xyz(coord, f"coord")

            # prepare scan setteing
            charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            nel = struc_inspector.get_value(f"{strucid}.struc_info.general.nel")

            if nel % 2 == 0:
                uhf = np.inf  # try both 0 and 2 (additional reaction path)
            else:
                uhf = 1

            ## constrain setting
            at1, at2 = info_l[1], info_l[2]
            at_l = [at1, at2]
            d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 2)
            d = searching_setting.get("pulling_distance", 2)
            push_l = [(int(at1) + 1,  # at1
                       int(at2) + 1,  # at2
                       d12,  # strat distance
                       d * d12)]  # end distance

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
            # call xtb scan
            status = CallxTB(scan_wd, {"method": xtb_setting.get("method", "gfn2"),
                                       "coordfn": "coord",
                                       "optlevel": xtb_setting.get("optlevel", None),
                                       "charge": charge,
                                       "uhf": uhf if not np.isinf(uhf) else 0,
                                       "alpb": xtb_setting.get("alpb", "h2o"),
                                       "etemp": xtb_setting.get("etemp", "4000")
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
                                         "atomic_connection": (constrain_l, ele)
                                         }  # endxyz info
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
                                           "etemp": xtb_setting.get("etemp", "4000")
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
                                             "atomic_connection": (constrain_l, ele)
                                             }  # endxyz info
                                            ]
                count += 1
        return output_d
    elif scanmode == 1:  # use molbar FF
        from aRST.script.do_scan import get_endxyz_from_mb
        from aRST.script.callsys import CallxTB
        for dgap in sorted(BO.keys())[:ndisso]:
            # prepare scan_wd
            scan_wd = os.path.join(wd, "scan", f"{reactionpath}_{count}")
            create_folder(scan_wd)

            # prepare scan_xyz
            info_l = BO[dgap][0]
            strucid = info_l[0]
            coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            nel = struc_inspector.get_value(f"{strucid}.struc_info.general.nel")
            charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            write_xyz(coord, f"coord")

            # get atomic connection
            atom_connection, ele = struc_inspector.get_value(f"{strucid}.struc_info.atomic_connection")
            for item in atom_connection:
                if item["atoms"] == info_l[1:3] or item["atoms"] == info_l[3:1]:
                    item["value"] *= 1.5

            coordinates = [list(coord.loc[line, 1:4]) for line in range(coord.shape[0])]
            mbxyz = get_endxyz_from_mb(atom_connection, coordinates, ele)
            write_xyz(mbxyz, "coord0")

            status = CallxTB(scan_wd, {"method": xtb_setting.get("method", "gfn2"),
                                       "coordfn": "coord0",
                                       "optlevel": xtb_setting.get("optlevel", None),
                                       "charge": charge,
                                       "uhf": 1 if nel % 2 != 0 else 2,
                                       "alpb": xtb_setting.get("alpb", "h2o"),
                                       "etemp": xtb_setting.get("etemp", "4000")
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
                                         "atomic_connection": (constrain_l, ele)
                                         }  # endxyz info
                                        ]
                count += 1
    else:
        raise Exception("unknown keyword for scan_mode, accecpted values include (0,1).")

    # # tmp for Mechnaophre min_path!!!!
    # if self.kwargs.get("terminalat") is not None:
    #     terminalat_setteing = self.kwargs.get("terminalat")
    #     strucid = "0001"
    #     at1 = terminalat_setteing.get("at1", None)
    #     at2 = terminalat_setteing.get("at2", None)
    #
    #     reaction_id = str(current_nreaction.value + 1).zfill(4)
    #     current_nreaction.value += 1
    #
    #     self.OUT_Reaction.set_value(f"{reaction_id}", {})  # occupy a position
    #     scan_wd = os.path.join(self.wd, "scan", f"{reaction_id}")
    #     create_folder(scan_wd)
    #
    #     coord = self.struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
    #     write_xyz(coord, "coord")
    #     d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 3)
    #     charge = self.struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
    #     uhf = self.struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1
    #
    #     push_l = [(str(int(at1) + 1),  # at1
    #                str(int(at2) + 1),  # at2
    #                str(d12),  # d_12
    #                str(d * d12))  # 2 * d_12
    #               ]
    #     scan_setting = {
    #         "constrain": {
    #             "k": xtb_setting.get("constrain_strength", str(1)),
    #             "distance": None
    #         },
    #         "scan": {
    #             "mode": None,
    #             "distance": push_l if push_l else None,
    #             "step": xtb_setting.get("scan_step", str(20))
    #         }
    #     }
    #     struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
    #     rc = struc_callxtb.call_scan(scan_setting)
    #     if rc == 0:
    #         str_pathnr = reactionpath + str(self.tmpcount)
    #         self.tmpcount += 1
    #         self.OUT_Reaction.set_value(f"{reaction_id}.reaction_path", [str_pathnr])
    #         print(f"resigning reaction id {reaction_id} to reaction path {str_pathnr}")
    #
    #         endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
    #         chargefile = os.path.join(scan_wd, "charges")
    #         nel = get_nelectron(os.path.join(scan_wd, "xtb.out"))
    #         products_l = self.followup_processsing_of_endxyz(endxyz, uhf + 1, charge, nel,
    #                                                          chargefile, current_nstruc)
    #         output_d[reaction_id] = [[strucid], products_l]
    #     return output_d
