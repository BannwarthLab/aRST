# Reaction simulation for intramolecular association
import json
from os.path import join
import shutil
import numpy as np
from aRST.script.toolbox import create_folder, write_xyz, read_coord, read_bo, read_atcharge, get_topo
from aRST.script.spr import SPR, SDP
from aRST.script.geom import Coord, alignment
from aRST.script.callsys import CallxTB
from aRST.script.do_scan import build_push_l
from aRST.script.do_scan import get_endxyz_from_mb
import aRST.script.initialization as storage
from molbar.barcode import get_molbar_from_coordinates
from aRST.script.head import print_react_l
import copy


def get_interA_react_l(reactants0):
    generalsetting = storage.setting.searching_setting.general
    intersetting = storage.setting.searching_setting.inter
    react_l = SPR(reactants0, generalsetting.spr_mode).inter(
        equivalent_afogap_tol=generalsetting.equivalent_afogap_tol,
        exclude_elements=generalsetting.exclude_elements,
        exclude_pairs=generalsetting.exclude_pairs
    )
    cut = 2 if intersetting.concerted_pair else 1
    for item in react_l:
        pairs = item["atl"][:cut]
        item["atl"] = pairs

    # print top 10 react_l
    print_react_l(reactants0, react_l)
    # write full react_l in file
    with open(join(generalsetting.wd, "react_l"), "w") as f:
        for item in react_l:
            f.write(json.dumps(item) + "\n")

    atom_pairs = intersetting.atom_pair
    if atom_pairs:
        # overwrite react_l with user's input
        # predifined reaction, used in one cycle
        at1, at2 = atom_pairs[0]
        react_l = [{"reactants": [reactants0[0], reactants0[1]],
                    "atl": [[at1, at2]],
                    "gap": 0,
                    "force_combination": True}]
        if len(atom_pairs) > 1:

            intersetting.concerted_pair = True
            for atom_pair in atom_pairs[1:]:
                at_combi_list = react_l[0]["atl"]
                at3, at4 = atom_pair
                at_combi_list.append([at3, at4])
    return react_l


def generate_reactive_complex(entry):
    intersetting = storage.setting.searching_setting.inter
    # prepare scan_xyz
    struc1, struc2 = entry["reactants"]
    s1 = storage.allstruc.get(struc1)
    coordinate1 = s1.general.coordinate
    ele_l1 = s1.general.ele_l
    nelA = s1.general.nel

    s2 = storage.allstruc.get(struc2)
    coordinate2 = s2.general.coordinate
    ele_l2 = s2.general.ele_l
    nelB = s2.general.nel
    nat_coord1 = len(ele_l1)

    # generate reactiv complex
    at_l = entry["atl"]
    coordinate_merge, ele_merged = alignment(coordinate1, ele_l1, coordinate2, ele_l2, at_l,
                                             extra_setting={"maxovlap": intersetting.maxovlap,
                                                            "aligning_distance": intersetting.aligning_distance})
    # prepare scan setteing
    charge = sum([k.general.charge for k in [s1, s2]])

    force_combination = False if "force_combination" not in entry.keys() else entry["force_combination"]
    push_l = build_push_l(coordinate_merge, ele_merged, at_l, nat_coord1, force_combination=force_combination
                          )
    return struc1, struc2, at_l, coordinate_merge, ele_merged, nat_coord1, charge, nelA + nelB, push_l


def check_scan_uhf_interA(entry):
    struc1, struc2 = entry["reactants"]
    nelA = storage.allstruc.get(struc1).general.nel
    nelB = storage.allstruc.get(struc2).general.nel
    if nelA % 2 == 0 and nelB % 2 == 0:
        return [0]
    elif (nelA + nelB) % 2 != 0:
        return [1]
    else:
        return [0, 2]


def interxtbscan(reactants0, entry, reactionpath, count=0):
    generalsetting = storage.setting.searching_setting.general
    intersetting = storage.setting.searching_setting.inter
    xtbsetting = storage.setting.xtb_setting
    scanwd = generalsetting.scanwd
    output_d = {}
    # preset output
    uhf_l = check_scan_uhf_interA(entry)
    struc1, struc2, at_l, \
    coordinate_merge, ele_merged, nat_coord1, \
    charge, totnel, push_l = generate_reactive_complex(entry)

    # preparing scan inp
    constrain_setting = {
        "constrain": {
            "force_constant": xtbsetting.k
        },
        "scan": {
            "mode": "concerted" if len(at_l) == 2 else None,
            "distance": push_l if push_l else None,
            "step": xtbsetting.step
        }
    }
    if intersetting.loose_bond:
        # searching weakest bond
        allbonds = SDP([struc1, struc2],
                       exclude_elements=generalsetting.exclude_elements,
                       exclude_pairs=generalsetting.exclude_pairs
                       ).allbonds
        allbonds = sorted(allbonds, key=lambda x: x["mbo"])
        item = allbonds[0]
        if item["mbo"] <= 0.9:
            if item["reactants"] == struc1:
                at1, at2 = item["atl"]
                at1, at2 = int(at1), int(at2)
            else:
                at1, at2 = item["atl"]
                at1, at2 = int(at1) + nat_coord1, int(at2) + nat_coord1
            dd = np.round(Coord(coordinate_merge, ele_merged).get_ats_distance(int(at1), int(at2)), 3)
            constrain_setting["constrain"]["loose_bond"] = at1 + 1, at2 + 1, dd

    unreacted = sorted((x for x in reactants0 if x != struc1 and x != struc2), key=int, reverse=True)
    for u in uhf_l:
        print(f"Simulating interA path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...", end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")
        create_folder(scan_wd)
        write_xyz(coordinate_merge, ele_merged, f"coord")

        with open(join(scan_wd, "atl"), "w") as f:
            f.write(json.dumps(entry))

        # call xtb scan
        status = CallxTB(scan_wd,
                         custom_settings={"charge": charge, "uhf": u,
                                          "etemp": xtbsetting.scan_etemp}).scan(constrain_setting)

        if status == 0:
            endxyz, ele_l_end = read_coord(join(scan_wd, "xtbopt.xyz"))
            bo = read_bo(join(scan_wd, "wbo"))
            constrain_l = get_topo(endxyz, ele_l_end, bo)
            atcharge = read_atcharge(join(scan_wd, "charges"))
            if generalsetting.sep_mode != 0:
                # do a unconstarined opt
                optwd = join(scan_wd, "opt")
                create_folder(optwd)
                write_xyz(endxyz, ele_l_end, "coord")
                status2 = CallxTB(optwd, custom_settings={"charge": charge, "uhf": u}).opt()
                if status2 == 0:
                    endxyz, ele_l_end = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l = get_topo(endxyz, ele_l_end, bo)
                    atcharge = read_atcharge(join(optwd, "charges"))
            # md optional
            if generalsetting.md:
                mdwd = join(scan_wd, "md")
                bo, constrain_l, atcharge, endxyz, ele_l_end = _short_md_for_product_struc_guess(mdwd, endxyz,
                                                                                                 ele_l_end, push_l,
                                                                                                 charge, u)

            endxyz_mb = get_molbar_from_coordinates(endxyz, ele_l_end)
            output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc1, struc2],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l, "simuhf": u,
                                                    "scanwd": f"{reactionpath}_{count}"},  # reactants info
                                                   {"atcharge": atcharge,
                                                    "charge": charge,
                                                    "nel": totnel,
                                                    "uhf": u,
                                                    "coordinate": endxyz,
                                                    "ele_l": ele_l_end,
                                                    "bo": bo,
                                                    "atomic_connection": constrain_l,
                                                    "mb":endxyz_mb},  # endxyz info
                                                   {"type": "interA"}  # reaction info
                                                   ]
            print(f"done.")
        else:
            print(f"failed.")
        count += 1
    return output_d, count


def intermbff(reactants0, entry, reactionpath, count=0):
    generalsetting = storage.setting.searching_setting.general
    intersetting = storage.setting.searching_setting.inter
    scanwd = generalsetting.scanwd
    output_d = {}
    uhf_l = check_scan_uhf_interA(entry)

    struc1, struc2, at_l, \
    coordinate_merge, ele_merged, nat_coord1, \
    charge, totnel, push_l = generate_reactive_complex(entry)

    # get atomic connection
    atom_connection1 = copy.deepcopy(storage.allstruc.get(struc1).general.atomic_connection)
    atom_connection2 = copy.deepcopy(storage.allstruc.get(struc2).general.atomic_connection)
    add_atom_connection2 = atom_connection2
    for item in add_atom_connection2:
        for at in range(len(item["atoms"])):
            item["atoms"][at] += nat_coord1

    atom_connection = atom_connection1 + add_atom_connection2
    for item in push_l:
        at1, at2 = int(item[0]), int(item[1])
        value = np.round(float(item[3]), 2)
        atom_connection.append({"atoms": [at1, at2], "value": value})

    unreacted = sorted((x for x in reactants0 if x != struc1 and x != struc2), key=int, reverse=True)
    mbxyz0, ele_l0 = get_endxyz_from_mb(atom_connection, coordinate_merge, ele_merged, set_angles=False)
    mbxyz, mbele_l = get_endxyz_from_mb(atom_connection, mbxyz0, ele_l0, set_angles=True)
    for u in uhf_l:
        print(f"Simulating interA path {reactionpath}_{count}: R={reactants0}, atl={at_l}, uhf_sim={u} ...", end="")
        scan_wd = join(scanwd, f"{reactionpath}_{count}")
        create_folder(scan_wd)
        write_xyz(coordinate_merge, ele_merged, f"coord")

        with open(join(scan_wd, "atl"), "w") as f:
            f.write(json.dumps(entry))

        with open(join(scan_wd, "newbond"), "w") as f:
            for item in atom_connection:
                f.write(json.dumps(item) + "\n")

        # mbxyz=get_endxyz_from_mb(atom_connection, coordinates, ele, set_angles=True)

        write_xyz(mbxyz0, ele_l0, "coord_mb0")
        write_xyz(mbxyz, mbele_l, "coord_mb")
        # constrained opt in gfn2
        status = CallxTB(scan_wd,
                         custom_settings={"charge": charge, "uhf": u,
                                          "k": storage.setting.xtb_setting.k,
                                          "coordfn": "coord_mb",
                                          "etemp": storage.setting.xtb_setting.scan_etemp
                                          }).constrained_opt(
            newbond=atom_connection[-len(push_l):])
        if status == 0:
            endxyz, ele_l_end = read_coord(join(scan_wd, "xtbopt.xyz"))
            bo = read_bo(join(scan_wd, "wbo"))
            constrain_l = get_topo(endxyz, ele_l_end, bo)
            atcharge = read_atcharge(join(scan_wd, "charges"))

            if generalsetting.sep_mode != 0:
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
            if generalsetting.md:
                mdwd = join(scan_wd, "md")
                bo, constrain_l, atcharge, endxyz, ele_l_end = _short_md_for_product_struc_guess(mdwd, endxyz,
                                                                                                 ele_l_end, push_l,
                                                                                                 charge, u)

            endxyz_mb = get_molbar_from_coordinates(endxyz, ele_l_end)
            output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc1, struc2],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l,
                                                    "simuhf": u,
                                                    "scanwd": f"{reactionpath}_{count}"},  # reactants info
                                                   {"atcharge": atcharge,
                                                    "charge": charge,
                                                    "nel": totnel,
                                                    "uhf": u,
                                                    "coordinate": endxyz,
                                                    "ele_l": ele_l_end,
                                                    "bo": bo,
                                                    "atomic_connection": constrain_l,
                                                    "mb":endxyz_mb},  # endxyz info
                                                   {"type": "interA"}  # reaction info
                                                   ]
            print(f"done.")
        else:
            print(f"failed.")
        count += 1
    return output_d, count


def _short_md_for_product_struc_guess(mdwd, coordinate, ele_l, constrained_bonds, charge, uhf):
    create_folder(mdwd)
    write_xyz(coordinate, ele_l, "coord")
    mdtemp = storage.setting.xtb_setting.scan_etemp
    mdsetting=storage.setting.searching_setting.md
    xtbmdparse = CallxTB(mdwd, custom_settings={"charge": charge, "uhf": uhf,
                                                "mdetemp": mdtemp,
                                                "mdtime": mdsetting.time,
                                                "mdstep": mdsetting.step,
                                                "mdshake": mdsetting.shake,
                                                "mddump": mdsetting.dump,
                                                "etemp": storage.setting.xtb_setting.scan_etemp
                                                # sp after md
                                                })
    status = xtbmdparse.short_md(constrained_bonds=constrained_bonds)
    if status != 0:
        raise ValueError(f"md has error, check {mdwd}")
    # get energy min struc in the md simulation
    # coord can't be single at in this case
    coordinate, ele_l, emd = xtbmdparse.get_Emin_geom_from_mdtrj(nat=len(ele_l))
    write_xyz(coordinate, ele_l, "emincoord")
    # doing a sp in the same folder for sp information
    status = xtbmdparse.sp()
    if status != 0:
        raise ValueError(f"sp after md has error, check {mdwd}")
    bo = read_bo(join(mdwd, "wbo"))
    constrain_l = get_topo(coordinate, ele_l, bo)
    atcharge = read_atcharge(join(mdwd, "charges"))
    return bo, constrain_l, atcharge, coordinate, ele_l
