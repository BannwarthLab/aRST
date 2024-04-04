# Reaction simulation for intramolecular association
import os

import numpy as np
import pandas as pd
from toolbox import create_folder, write_xyz
from spr import SPR_dissociation

from geom import Coord

def simulation_intraD(reactants, reactionpath, struc_inspector, searching_setting, xtb_setting, ndisso, wd):
    """
    :param reactants: reactants list for this step
    :param reactionpath: str index of rp
    :param struc_inspector: struc info d
    :param searching_setting: typed
    :param xtb_setting: typed
    :param ndisso: number of searching range
    :param wd: working dictory
    reactants_l: reactants list really react
    :return: output_d[count] = (reactants_l, at_l, gap, endxyz, chargefile, nel)
    """

    output_d = {}



    ###### rangking reactivity ######
    BO = SPR_dissociation(reactants,
                          struc_inspector,
                          allow_hydrogen=searching_setting.get("hydrogen"),
                          type=searching_setting.get("dissotype", "wboe"),
                          allow_halogen=searching_setting.get("halogen"))

    if not BO:
        return output_d
    print("disassoication possibility")
    for key in sorted(BO.keys())[:10]:
        print(key, BO[key])

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





    ###### perform simulation ######
    if xtb_setting.get('scanmode', None) is None or xtb_setting.get('scanmode', None) == "xtb":
        from callsys import CallxTB
        count = 1
        # standard scan process
        for dgap in sorted(BO.keys())[:ndisso]:
            # prepare scan_wd
            scan_wd = os.path.join(wd, "scan", f"intraD_{reactionpath}{count}")
            create_folder(scan_wd)

            # prepare scan_xyz
            info_l = BO[dgap][0]
            strucid = info_l[0].zfill(4)
            coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            write_xyz(coord, f"coord")

            # prepare scan setteing
            ## global features
            charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            uhf = struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1

            ## constrain setting
            at1,at2 = info_l[1],info_l[2]
            at_l = [at1,at2]

            d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 2)
            d = searching_setting.get("ddissociation", 2)
            push_l = [(str(int(at1) + 1),  # at1
                       str(int(at2) + 1),  # at2
                       str(d12),  # strat distance
                       str(d * d12))]  #  end distance
            scan_setting = {
                "constrain":
                    {"k": xtb_setting.get("constrain_strength", str(1)),
                     "distance": None},
                "scan":
                    {"mode": "concerted" if len(at_l) == 2 else None,
                     "distance": push_l if push_l else None,
                     "step": xtb_setting.get("scan_step", str(20))}
            }

            # call xtb scan
            level = xtb_setting.get('xtblevel', 'gfn2')
            alpb = xtb_setting.get('alpb', None)
            etemp = xtb_setting.get("etemp", str(300))
            struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
            status = struc_callxtb.call_scan(scan_setting)


            if status ==0:
                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                chargefile = os.path.join(scan_wd, "charges")
                nel = struc_inspector.get_value(f"{strucid}.struc_info.general.nel")

                output_d[count] = ([strucid], at_l, dgap, endxyz, chargefile, nel)
            count += 1
        return output_d
    else: # use mb ff
        from do_scan import get_endxyz_from_mb, evaluate_bo_after_mbff
        for dgap in sorted(BO.keys())[:ndisso]:
            info_l = BO[dgap][0]
            strucid = info_l[0].zfill(4)
            at1,at2 = info_l[1], info_l[2]
            atom_connection, ele = struc_inspector.get_value(f"{strucid}.struc_info.atomic_connection")

            scan_wd = os.path.join(wd, "scan", f"intraD_{reactionpath}{count}")
            create_folder(scan_wd)

            coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            write_xyz(coord, f"coord")
            with open(f"constrainl", "w") as f:
                f.write(str(at1)+"," +str(at2) +  "\n")

            coordinates = [coord.loc[line, 1:4].tolist() for line in range(coord.shape[0])]

            atom_connection.append({"atoms": [at1+1, at2+1],
                                    "value": d*(Coord(coord).get_rr_distance(at1) +
                                             Coord(coord).get_rr_distance(at2))})

            endxyz = get_endxyz_from_mb(atom_connection, coordinates, ele)
            write_xyz(endxyz, f"endxyz")

            charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            uhf = struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1
            nel = struc_inspector.get_value(f"{strucid}.struc_info.general.nel")
            status = evaluate_bo_after_mbff(scan_wd, level, endxyz, charge, uhf, alpb, etemp, None)
            if status:
                chargefile = os.path.join(scan_wd, "charges")
                output_d[count] = ([strucid], [at1,at2], dgap, endxyz, chargefile, nel)
            count+=1
        return output_d
