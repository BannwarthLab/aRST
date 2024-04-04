# Reaction simulation for intramolecular association
import os

import numpy as np
import pandas as pd
from toolbox import create_folder, write_xyz
from spr import SPR_HL, SPR_dissociation

from geom import Coord, alignment




def simulation_interA(reactants, reactionpath, struc_inspector, searching_setting, xtb_setting, ninter, wd):
    """
    :param reactants: reactants list for this step
    :param reactionpath: str index of rp
    :param struc_inspector: struc info d
    :param searching_setting: typed
    :param xtb_setting: typed
    :param ninter: number of searching range
    :param wd: working dictory
    reactants_l: reactants list really react
    :return: output_d[count] = (reactants_l, at_l, gap, endxyz, chargefile, nel)
    """
    output_d = {}

    ###### rangking reactivity ######
    react_os, react_ss, HOAO_a, HOAO_b, LUAO_a, LUAO_b = SPR_HL(reactants,
                                                                struc_inspector,
                                                                allow_hydrogen=searching_setting.get("hydrogen"),
                                                                search_range=2 if searching_setting.get(
                                                                    "samestruc") else 1)
    BO = SPR_dissociation(reactants,
                          struc_inspector,
                          allow_hydrogen=searching_setting.get("hydrogen"),
                          type=searching_setting.get("dissotype", "wboe"),
                          allow_halogen=searching_setting.get("halogen"))

    intera = {}
    for react_d in [react_os, react_ss]:
        for current_key, current_item in react_d.items():
            # if current_key>0: # ? discard hl<0
            #  5.93 {('1', '1'): [['0', '0']]}
            for struc_combi, at_l in current_item.items():
                struc1, struc2 = struc_combi
                intera[current_key] = struc1, struc2, at_l

    # react_d:  key= 5.93, value = {('1', '1'), [['0', '0']]}
    # intera:   key= 5.93, value = ('1','1',[['0', '0']])
    print("inter assoication possibility")
    for key in sorted(intera.keys())[:10]:
        print(key, intera[key])



    from dictionaryInspector import remove_idential_atl
    intera = remove_idential_atl(intera)
    intera = {key: intera[key] for key in sorted(intera.keys())[:10]}

    if os.path.isfile(os.path.join(wd,"atl")):
        os.remove(os.path.join(wd,"atl"))
    with open(os.path.join(wd,"atl"),"a+") as f:
        f.write(f"nr,gap,at_l,bo,bl,bl/covalentradii\n")



    if not intera:
        return output_d

    def constrain_weakest_bond(struc1,struc2, BO, bolim):
        def get_weakest_bo(struc):
            strucid = struc.rstrip("0")
            tempBO = {key: [value[0][1], value[0][2]] for key, value in BO.items() if value[0][0] == strucid}
            if not tempBO:
                return None, None, None
            lowest_bo = min(tempBO.keys())
            at_l = tempBO[lowest_bo]
            at1, at2 = at_l[0], at_l[1]
            return lowest_bo, at1, at2

        tmpr, tmpat1, tmpat2 = None, None, None
        lowest_bo1, tmp1at1, tmp1at2 = get_weakest_bo(struc1)
        tmpr1 = np.round(Coord(coord1).get_ats_distance(int(tmp1at1), int(tmp1at2)), 3) \
            if tmp1at1 is not None and tmp1at2 is not None else None

        lowest_bo2, tmp2at1, tmp2at2 = get_weakest_bo(struc2)
        if tmp2at1 is not None and tmp2at2 is not None:
            tmpr2 = np.round(Coord(coord2).get_ats_distance(int(tmp2at1), int(tmp2at2)), 3)
            tmp2at1, tmp2at2 = str(int(tmp2at1) + nat_coord1), str(int(tmp2at2) + nat_coord1)
        else:
            tmpr2, tmp2at1, tmp2at2 = None, None, None

        lowest_bo1 = np.inf if lowest_bo1 is None else lowest_bo1
        lowest_bo2 = np.inf if lowest_bo2 is None else lowest_bo2

        if lowest_bo2 < bolim and lowest_bo2 < lowest_bo1:
            tmpr = tmpr2
            tmpat1, tmpat2 = tmp2at1, tmp2at2
        elif lowest_bo1 < bolim and lowest_bo1 < lowest_bo2:
            tmpr = tmpr1
            tmpat1, tmpat2 = tmp1at1, tmp1at2
        return tmpat1, tmpat2,tmpr


    ###### perform simulation ######
    from do_scan import build_push_l
    # global xtb setting
    level = xtb_setting.get('xtblevel', 'gfn2')
    alpb = xtb_setting.get('alpb', None)
    etemp = xtb_setting.get("etemp", str(300))

    if xtb_setting.get('scanmode', None) is None or xtb_setting.get('scanmode', None) == "xtb":
        from callsys import CallxTB

        # standard scan process

        intera = {key: intera[key] for key in sorted(intera.keys())[:ninter]}

        count = 1
        for gap, items in intera.items():

            # prepare scan_wd
            scan_wd = os.path.join(wd, "scan", f"interA_{reactionpath}{count}")
            create_folder(scan_wd)

            # prepare scan_xyz
            struc1, struc2 = intera[gap][0], intera[gap][1]
            struc1, struc2 = struc1.zfill(4), struc2.zfill(4)

            coord1 = struc_inspector.get_value(f"{struc1}.struc_info.general.coord")
            coord2 = struc_inspector.get_value(f"{struc2}.struc_info.general.coord")
            nat_coord1 = int(coord1.shape[0])

            ## alignment according to first atomic pair
            at_l = intera[gap][2]
            at1, at2 = int(at_l[0][0]), int(at_l[0][1])
            alignedcoord = alignment(coord1, coord2, at1, at2, **xtb_setting)
            write_xyz(alignedcoord, f"coord")

            # prepare scan setteing
            ## global features
            charge = sum([struc_inspector.get_value(f"{k}.struc_info.general.charge") for k in [struc1, struc2]])
            uhf = sum([struc_inspector.get_value(f"{k}.struc_info.general.multiplicity") for k in [struc1, struc2]]) - 2

            ## constrain setting
            push_l = build_push_l(alignedcoord, at_l,nat_coord1)

            ## loose the weakest bond if bo < 0.9
            weakestbo_at1, weakestbo_at2, weakestbo_r = constrain_weakest_bond(struc1,struc2,BO,0.9)
            scan_setting = {
                "constrain":
                    {"k": xtb_setting.get("constrain_strength", str(1)),
                     "distance": [str(weakestbo_at1 + 1),
                                  str(weakestbo_at2 + 1),
                                  str(weakestbo_r * 1.5)] if weakestbo_r is not None else None
                     },
                "scan":
                    {"mode": "concerted" if len(at_l) == 2 else None,
                     "distance": push_l if push_l else None,
                     "step": xtb_setting.get("scan_step", str(20))}
            }

            # call xtb scan
            struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
            status = struc_callxtb.call_scan(scan_setting)

            if status ==0:
                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                chargefile = os.path.join(scan_wd, "charges")
                nel = sum([struc_inspector.get_value(f"{k}.struc_info.general.nel") for k in [struc1,struc2]])
                output_d[count] = ([struc1,struc2], at_l, gap, endxyz, chargefile, nel)

            count += 1
        return output_d
    else:  # use mb ff to quick screen all 10 options
        from do_scan import get_endxyz_from_mb, evaluate_bo_after_mbff,get_bo,get_bl

        # prepare scan_wd (all in one)
        scan_wd = os.path.join(wd, "scan", f"interA_{reactionpath}")
        create_folder(scan_wd)
        count = 1
        tmp=1
        tmp_value=None
        for gap, items in intera.items():
            os.chdir(scan_wd)
            # prepare scan_xyz
            struc1, struc2 = intera[gap][0], intera[gap][1]
            struc1, struc2 = struc1.zfill(4), struc2.zfill(4)

            coord1 = struc_inspector.get_value(f"{struc1}.struc_info.general.coord")
            coord2 = struc_inspector.get_value(f"{struc2}.struc_info.general.coord")
            nat_coord1 = int(coord1.shape[0])

            ## alignment according to first atomic pair
            at_l = intera[gap][2]
            at1, at2 = int(at_l[0][0]), int(at_l[0][1])
            alignedcoord = alignment(coord1, coord2, at1, at2, **xtb_setting)
            write_xyz(alignedcoord, f"coord_{count}")

            # prepare scan setteing
            ## global features
            charge = sum([struc_inspector.get_value(f"{k}.struc_info.general.charge") for k in [struc1, struc2]])
            uhf = sum([struc_inspector.get_value(f"{k}.struc_info.general.multiplicity") for k in [struc1, struc2]]) - 2
            nel = sum([struc_inspector.get_value(f"{k}.struc_info.general.nel") for k in [struc1, struc2]])

            ## constrain setting
            push_l = build_push_l(alignedcoord, at_l, nat_coord1)


            atom_connection1, ele1 = struc_inspector.get_value(f"{struc1.zfill(4)}.struc_info.atomic_connection")
            atom_connection2, ele2 = struc_inspector.get_value(f"{struc2.zfill(4)}.struc_info.atomic_connection",
                                                               nat_coord1)
            # atom_connection = atom_connection1 + atom_connection2
            atom_connection = []
            ele = ele1 + ele2

            new_connection = []
            for item in push_l:
                at1,at2 = int(item[0]),int(item[1])
                value = float(item[3])
                tmp_value=value
                atom_connection.append({"atoms": [at1, at2],
                                        "value": value})

                new_connection.append([at1,at2])
                with open(f"constrainl_{count}", "a+") as f:
                    for atl in atom_connection:
                        f.write(f"{atl}\n")

            ## loose the weakest bond if bo < 0.9
            weakestbo_at1, weakestbo_at2, weakestbo_r = constrain_weakest_bond(struc1, struc2, BO, 0.9)
            if weakestbo_r is not None:
                atom_connection.append({"atoms": [int(weakestbo_at1)+1, int(weakestbo_at2)+1],
                                        "value": float(weakestbo_r)})


            coordinates = [list(alignedcoord.loc[line, 1:4]) for line in range(alignedcoord.shape[0])]
            endxyz = get_endxyz_from_mb(atom_connection, coordinates, ele)
            write_xyz(endxyz, f"endxyz_{count}")

            # evaluate new created bond with wbo
            cwd = os.path.join(scan_wd, str(count))
            create_folder(cwd)

            status = evaluate_bo_after_mbff(cwd, level, endxyz, charge, uhf, alpb, etemp, new_connection)
            endxyz = pd.read_table(os.path.join(cwd, "xtbopt.xyz"), header=None, skiprows=2, sep='\s+')

            bo = get_bo(cwd,new_connection)
            bl = get_bl(endxyz,new_connection)

            if status:
                chargefile = os.path.join(cwd, "charges")
                output_d[count] = ([struc1, struc2], at_l, gap, endxyz, chargefile, nel)

            with open(os.path.join(wd, "atl"), "a+") as f:
                f.write(f"{tmp},{gap},{intera[gap]},{bo},{bl},{np.round(bl*0.95/tmp_value,2)}\n")

            count += 1
            tmp += 1


    return output_d


def force_interA(reactants, reactionpath, struc_inspector, forceconstrain_setting, xtb_setting, ninter, wd):
    """
    :param reactants: reactants list for this step
    :param reactionpath: str index of rp
    :param struc_inspector: struc info d
    :param forceconstrain_setting: typed
    :param xtb_setting: typed
    :param ninter: number of searching range
    :param wd: working dictory
    reactants_l: reactants list really react
    :return: output_d[count] = (reactants_l, at_l, gap, endxyz, chargefile, nel)
    """
    output_d = {}

    ###### perform simulation ######
    from do_scan import build_push_l
    # global xtb setting
    level = xtb_setting.get('xtblevel', 'gfn2')
    alpb = xtb_setting.get('alpb', None)
    etemp = xtb_setting.get("etemp", str(300))

    if xtb_setting.get('scanmode', None) is None or xtb_setting.get('scanmode', None) == "xtb":
        from callsys import CallxTB

        # standard scan process

        count = 1
        for gap, items in forceconstrain_setting.items():

            # prepare scan_wd
            scan_wd = os.path.join(wd, "scan", f"interA_{reactionpath}{count}")
            create_folder(scan_wd)

            # prepare scan_xyz
            struc1, struc2 = str(items["stru1"]),str(items["stru2"])
            struc1, struc2 = struc1.zfill(4), struc2.zfill(4)

            coord1 = struc_inspector.get_value(f"{struc1}.struc_info.general.coord")
            coord2 = struc_inspector.get_value(f"{struc2}.struc_info.general.coord")
            nat_coord1 = int(coord1.shape[0])

            ## alignment according to first atomic pair
            at1, at2 = items["at1"],items["at2"]
            alignedcoord = alignment(coord1, coord2, at1, at2, **xtb_setting)
            write_xyz(alignedcoord, f"coord")

            # prepare scan setteing
            ## global features
            charge = sum([struc_inspector.get_value(f"{k}.struc_info.general.charge") for k in [struc1, struc2]])
            uhf = sum([struc_inspector.get_value(f"{k}.struc_info.general.multiplicity") for k in [struc1, struc2]]) - 2

            ## constrain setting
            if "at4" in items:
                at3,at4 = items["at3"],items["at4"]
                at_l = [[f'{at1}', f'{at2}'],[f'{at3}', f'{at4}']]
            else:
                at_l = [[f'{at1}', f'{at2}']]
            push_l = build_push_l(alignedcoord, at_l, nat_coord1)

            ## loose the weakest bond if bo < 0.9
            scan_setting = {
                "constrain":
                    {"k": xtb_setting.get("constrain_strength", str(1)),
                     "distance":None
                     },
                "scan":
                    {"mode": "concerted" if len(at_l) == 2 else None,
                     "distance": push_l if push_l else None,
                     "step": xtb_setting.get("scan_step", str(20))}
            }

            # call xtb scan
            struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
            status = struc_callxtb.call_scan(scan_setting)

            if status == 0:
                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                chargefile = os.path.join(scan_wd, "charges")
                nel = sum([struc_inspector.get_value(f"{k}.struc_info.general.nel") for k in [struc1, struc2]])
                output_d[count] = ([struc1, struc2], at_l, gap, endxyz, chargefile, nel)

            count += 1
        return output_d
    # else:  # use mb ff to quick screen all 10 options
    #     from do_scan import get_endxyz_from_mb, evaluate_bo_after_mbff, get_bo, get_bl
    #
    #     # prepare scan_wd (all in one)
    #     scan_wd = os.path.join(wd, "scan", f"interA_{reactionpath}")
    #     create_folder(scan_wd)
    #     count = 1
    #     tmp = 1
    #     tmp_value = None
    #     for gap, items in intera.items():
    #         os.chdir(scan_wd)
    #         # prepare scan_xyz
    #         struc1, struc2 = intera[gap][0], intera[gap][1]
    #         struc1, struc2 = struc1.zfill(4), struc2.zfill(4)
    #
    #         coord1 = struc_inspector.get_value(f"{struc1}.struc_info.general.coord")
    #         coord2 = struc_inspector.get_value(f"{struc2}.struc_info.general.coord")
    #         nat_coord1 = int(coord1.shape[0])
    #
    #         ## alignment according to first atomic pair
    #         at_l = intera[gap][2]
    #         at1, at2 = int(at_l[0][0]), int(at_l[0][1])
    #         alignedcoord = alignment(coord1, coord2, at1, at2, **xtb_setting)
    #         write_xyz(alignedcoord, f"coord_{count}")
    #
    #         # prepare scan setteing
    #         ## global features
    #         charge = sum([struc_inspector.get_value(f"{k}.struc_info.general.charge") for k in [struc1, struc2]])
    #         uhf = sum([struc_inspector.get_value(f"{k}.struc_info.general.multiplicity") for k in [struc1, struc2]]) - 2
    #         nel = sum([struc_inspector.get_value(f"{k}.struc_info.general.nel") for k in [struc1, struc2]])
    #
    #         ## constrain setting
    #         push_l = build_push_l(alignedcoord, at_l, nat_coord1)
    #
    #         atom_connection1, ele1 = struc_inspector.get_value(f"{struc1.zfill(4)}.struc_info.atomic_connection")
    #         atom_connection2, ele2 = struc_inspector.get_value(f"{struc2.zfill(4)}.struc_info.atomic_connection",
    #                                                            nat_coord1)
    #         # atom_connection = atom_connection1 + atom_connection2
    #         atom_connection = []
    #         ele = ele1 + ele2
    #
    #         new_connection = []
    #         for item in push_l:
    #             at1, at2 = int(item[0]), int(item[1])
    #             value = float(item[3])
    #             tmp_value = value
    #             atom_connection.append({"atoms": [at1, at2],
    #                                     "value": value})
    #
    #             new_connection.append([at1, at2])
    #             with open(f"constrainl_{count}", "a+") as f:
    #                 for atl in atom_connection:
    #                     f.write(f"{atl}\n")
    #
    #         ## loose the weakest bond if bo < 0.9
    #         weakestbo_at1, weakestbo_at2, weakestbo_r = constrain_weakest_bond(struc1, struc2, BO, 0.9)
    #         if weakestbo_r is not None:
    #             atom_connection.append({"atoms": [int(weakestbo_at1) + 1, int(weakestbo_at2) + 1],
    #                                     "value": float(weakestbo_r)})
    #
    #         coordinates = [list(alignedcoord.loc[line, 1:4]) for line in range(alignedcoord.shape[0])]
    #         endxyz = get_endxyz_from_mb(atom_connection, coordinates, ele)
    #         write_xyz(endxyz, f"endxyz_{count}")
    #
    #         # evaluate new created bond with wbo
    #         cwd = os.path.join(scan_wd, str(count))
    #         create_folder(cwd)
    #
    #         status = evaluate_bo_after_mbff(cwd, level, endxyz, charge, uhf, alpb, etemp, new_connection)
    #         endxyz = pd.read_table(os.path.join(cwd, "xtbopt.xyz"), header=None, skiprows=2, sep='\s+')
    #
    #         bo = get_bo(cwd, new_connection)
    #         bl = get_bl(endxyz, new_connection)
    #
    #         if status:
    #             chargefile = os.path.join(cwd, "charges")
    #             output_d[count] = ([struc1, struc2], at_l, gap, endxyz, chargefile, nel)
    #
    #         with open(os.path.join(wd, "atl"), "a+") as f:
    #             f.write(f"{tmp},{gap},{intera[gap]},{bo},{bl},{np.round(bl * 0.95 / tmp_value, 2)}\n")
    #
    #         count += 1
    #         tmp += 1

    return output_d
