# Reaction simulation for intramolecular association
import os
import shutil
import json
import numpy as np
import pandas as pd
from aRST.script.toolbox import create_folder, write_xyz, read_topomol, round_charge
from aRST.script.spr import SPR_dissociation
from aRST.script.dictionaryInspector import DictionaryInspector
from aRST.script.geom import Coord
from aRST.script.callsys import CallxTB


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
    ovwrite_at1 = searching_setting.get("at1", False)
    ovwrite_at2 = searching_setting.get("at2", False)


    #todo
    # react_l = SPR_dissociation(reactants,  # reactants id list
    #                            struc_rec,  # struc dict record
    #                            allow_hydrogen=searching_setting.get("hydrogen", True),
    #                            allow_halogen=searching_setting.get("halogen", True),
    #                            search_mode=searching_setting.get("mode", 1))
    react_l=[]
    if ovwrite_at1:
        if searching_setting.get("at3", False):
            ovwrite_at3 = searching_setting.get("at3", False)
            ovwrite_at4 = searching_setting.get("at4", False)
            react_l = [{"struc": str(reactants[0]),
                        "at_combi": [[str(ovwrite_at1), str(ovwrite_at2)], [str(ovwrite_at3), str(ovwrite_at4)]],
                        "dgap": None}]

        else:react_l = [{"struc": str(reactants[0]),
                        "at_combi": [[str(ovwrite_at1), str(ovwrite_at2)]],
                        "dgap": None}]



        if not react_l:
            raise Exception(f"no possible dissociation possibilities can be found in reactants list: {reactants}.")

    scanmode = searching_setting.get("scan_mode", 0)
    ndisso = searching_setting.get("number", 1)

    if scanmode == 0:
        # standard xtbscan

        for entry in react_l[:ndisso]:
            scan_wd = os.path.join(wd, "scan", f"{reactionpath}_{count}")
            create_folder(scan_wd)

            with open(os.path.join(scan_wd, "atl"), "w") as f:
                f.write(json.dumps(entry))

            # prepare scan_xyz
            struc = entry["struc"]
            coord = struc_inspector.get_value(f"{struc}.struc_info.general.coord")
            write_xyz(coord, f"coord")
            totat_l = entry["at_combi"]  # index starting with 1

            # prepare scan
            charge = struc_inspector.get_value(f"{struc}.struc_info.general.charge")
            nel = struc_inspector.get_value(f"{struc}.struc_info.general.nel")
            uhf=0
            # uhf = struc_inspector.get_value(f"{struc}.struc_info.general.uhf") # todo

            # preparing scan inp
            d = searching_setting.get("pulling_distance", 2)
            push_l = []
            for pair in totat_l[:2]:
                at1, at2 = pair
                d12 = np.round(Coord(coord).get_ats_distance(int(at1) - 1, int(at2) - 1), 2)
                dend = d * d12
                push_l.append([at1, at2, d12, dend])

            constrain_setting = {
                "constrain": {
                    "force_constant": xtb_setting.get("k", 1)
                },
                "scan": {
                    "mode": "concerted" if len(push_l) == 2 else None,
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
                                              "uhf": uhf if not np.isinf(uhf) else 0,
                                              "alpb": xtb_setting.get("alpb", "h2o"),
                                              "etemp": xtb_setting.get("struc_etemp", "4000")
                                              }).opt()
                    if status2 == 0:
                        endwd = optwd
                endxyz = pd.read_table(os.path.join(endwd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                if os.path.exists(os.path.join(endwd, 'xtbtopo.mol')):
                    constrain_l, ele = read_topomol(os.path.join(endwd, 'xtbtopo.mol'))
                else:
                    constrain_l, ele = [], []
                bo = pd.read_csv(os.path.join(endwd, "wbo"), header=None, sep='\s+') if os.path.getsize(
                    os.path.join(endwd, "wbo")) != 0 else None
                atcharge = pd.read_table(os.path.join(endwd, "charges"), header=None, sep='\s+')

                output_d[f"{count}"] = [{"reactants": [struc],
                                         "unreacted": [st for st in reactants
                                                       if st != struc]
                                         if len(reactants) > 2 else []},  # reactants info
                                        {"atomic_charge": atcharge,
                                         "charge": charge,
                                         "nel": nel,
                                         "uhf": uhf if not np.isinf(uhf) else 0,
                                         "coord": endxyz,
                                         "bo": bo,
                                         "atomic_connection": (constrain_l, ele)}  # endxyz info
                                        ]


    elif scanmode == 1:  # use molbar FF
        from aRST.script.do_scan import get_endxyz_from_mb
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
            atom_connection, ele = struc_inspector.get_value(f"{strucid}.struc_info.general.atomic_connection")
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
                                       "etemp": xtb_setting.get("scan_etemp", "4000")
                                       }).opt()
            if status == 0:
                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                if os.path.exists(os.path.join(scan_wd, 'xtbtopo.mol')):
                    constrain_l, ele = read_topomol(os.path.join(scan_wd, 'xtbtopo.mol'))
                else:
                    constrain_l, ele = [], []
                if os.path.getsize(os.path.join(scan_wd, "wbo")) != 0:
                    bo = pd.read_csv(os.path.join(scan_wd, "wbo"), header=None, sep='\s+')
                else:
                    bo = None
                atcharge = pd.read_table(os.path.join(scan_wd, "charges"), header=None, sep='\s+')
                output_d[f"{count}"] = [{"reactants": [strucid],
                                         "unreacted": [st for st in reactants if st != strucid]
                                         if len(reactants) > 1 else []},  # reactants info
                                        {"atomic_charge": atcharge,
                                         "charge": charge,
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
