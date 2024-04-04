# Reaction simulation for intramolecular association
import os
import pandas as pd
from toolbox import create_folder, write_xyz
from spr import SPR_HL, build_bo_d

from geom import Coord



def simulation_intraA(reactants, reactionpath, struc_inspector, searching_setting, xtb_setting, nintraa, wd):
    """
    :param reactants: reactants list for this step
    :param reactionpath: str index of rp
    :param struc_inspector: struc info d
    :param searching_setting: typed
    :param xtb_setting: typed
    :param nintraa: number of searching range
    :param wd: working dictory
    reactants_l: reactants list really react
    :return: output_d[count] = (reactants_l, at_l, gap, endxyz, chargefile, nel)
    """
    output_d = {}

    ###### rangking reactivity ######
    react_os, react_ss, HOAO_a, HOAO_b, LUAO_a, LUAO_b = SPR_HL(reactants,
                                                                struc_inspector,
                                                                allow_hydrogen=searching_setting.get("hydrogen"),
                                                                search_range=0)
    bo = build_bo_d(reactants,
                    struc_inspector,
                    allow_hydrogen=searching_setting.get("hydrogen"))  # bond connection info
    if not bo:
        return output_d

    intraa = {}
    for react_d in [react_os, react_ss]:
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
    print("intra assoication possibility")
    for key in sorted(intraa.keys())[:10]:
        print(key, intraa[key])

    intraa = {key: value for key, value in sorted(intraa.items())[:10]}


    if not intraa:
        return output_d


    ###### perform simulation ######
    from do_scan import build_push_l
    if xtb_setting.get('scanmode', None) is None or xtb_setting.get('scanmode', None) == "xtb":
        from callsys import CallxTB
        # global xtb setting
        level = xtb_setting.get('xtblevel', 'gfn2')
        alpb = xtb_setting.get('alpb', None)
        etemp = xtb_setting.get("etemp", str(300))


        # standard scan process
        intraa = {key: value for key, value in sorted(intraa.items())[:nintraa]}
        count=1
        for gap, items in intraa.items():

            # prepare scan_wd
            scan_wd = os.path.join(wd, "scan", f"intraA_{reactionpath}{count}")
            create_folder(scan_wd)

            # prepare scan_xyz
            strucid = intraa[gap][0].zfill(4)
            coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            write_xyz(coord, f"coord")

            # prepare scan setteing
            ## global features
            charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            uhf = struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1

            ## constrain setting
            at_l = intraa[gap][1]
            push_l = build_push_l(coord,at_l)
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
            struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
            status = struc_callxtb.call_scan(scan_setting)


            if status ==0:
                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                chargefile = os.path.join(scan_wd, "charges")
                nel = struc_inspector.get_value(f"{strucid}.struc_info.general.nel")

                output_d[count] = ([strucid],at_l,gap,endxyz,chargefile,nel)

            count += 1
        return output_d

    else: # use mb ff to quick screen all 10 options
        from do_scan import get_endxyz_from_mb,evaluate_bo_after_mbff

        scan_wd = os.path.join(wd, "scan", f"intraA_{reactionpath}")
        create_folder(scan_wd)

        count=1
        for gap, items in intraa.items():
            at_l = intraa[gap][1]
            strucid = intraa[gap][0].zfill(4)

            coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            write_xyz(coord, f"coord_{count}")

            if len(at_l) == 2:
                from do_scan import rerange_atl
                at1, at2, at3, at4 = int(at_l[0][0]), int(at_l[0][1]), \
                                     int(at_l[1][0]), int(at_l[1][1])
                at_l = rerange_atl(coord, at1, at2, at3, at4)


            atom_connection, ele = struc_inspector.get_value(f"{strucid}.struc_info.atomic_connection")

            with open(f"constrainl_{count}", "w") as f:
                for item in at_l:
                    f.write(str(item) + "\n")
            coordinates = [coord.loc[line, 1:4].tolist() for line in range(coord.shape[0])]

            for item in at_l:
                at1,at2= item
                atom_connection.append({"atoms": [at1+1, at2+1],
                                        "value": Coord(coord).get_rr_distance(at1) +
                                                Coord(coord).get_rr_distance(at2)})


            endxyz = get_endxyz_from_mb(atom_connection, coordinates, ele)
            write_xyz(endxyz, f"endxyz_{count}")

            # evaluate new created bond with wbo
            cwd = os.path.join(scan_wd,str(count))
            create_folder(cwd)
            charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            uhf = struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1
            nel = struc_inspector.get_value(f"{strucid}.struc_info.general.nel")
            status = evaluate_bo_after_mbff(cwd,level, endxyz, charge, uhf, alpb, etemp, at_l)

            if status:
                chargefile=os.path.join(cwd, "charges")
                output_d[count] = ([strucid], at_l, gap, endxyz, chargefile, nel)
                count+=1

        return output_d