import aRST.script.initialization as storage
from aRST.script.exploration.assign_reaction import SingleCycleSearch
from aRST.script.spr import SPR
from aRST.script.reactivity.print_stuff import *
from aRST.script.toolbox import create_folder, write_xyz
from aRST.script.geom import alignment


def aligning_two_structures(reactants=None, atl=None):
    if reactants is None:
        reactants = storage.allstruc.reactants_ini
    if len(reactants) != 2:
        raise ValueError("The number of input structures is not 2, please check")
    if atl is None:
        # check if define reactive ats
        atom_pairs = storage.setting.get_value("searching_setting.inter.atom_pair")
        if atom_pairs is not None:
            # overwrite react_l with user's input
            # predifined reaction, used in one cycle
            at1, at2 = atom_pairs[0]
            react_l = [{"struc_combi": [reactants[0], reactants[1]],
                        "at_combi": [[str(at1), str(at2)]],
                        "gap": 0,
                        "force_combination": True}]
            if len(atom_pairs) > 1:
                storage.setting.set_value("searching_setting.inter.concerted_pair", True)
                for atom_pair in atom_pairs[1:]:
                    at_combi_list = react_l[0]["at_combi"]
                    at3, at4 = atom_pair
                    at_combi_list.append([str(at3), str(at4)])
        else:
            print("aRST doesn't detect preassign reactive ats, switch to most reactive pair assigning by AFO gaps")
            SingleCycleSearch(reactants)
            spr = SPR(reactants, storage.setting.get_value("searching_setting.general.spr_mode", 0))
            react_l = spr.inter(
                equivalent_afogap_tol=storage.setting.get_value("searching_setting.general.equivalent_afogap_tol",
                                                                0.01),
                exclude_elements=storage.setting.get_value("searching_setting.general.exclude_elements"),
                exclude_pairs=storage.setting.get_value("searching_setting.general.exclude_pairs"))
            cut = 2 if storage.setting.get_value("searching_setting.inter.concerted_pair", False) else 1
            for item in react_l:
                pairs = item["at_combi"][:cut]
                item["at_combi"] = pairs
            entry = react_l[0]
            atl = entry["at_combi"]

    align_wd = join(storage.setting.wd0, f"alignment")
    create_folder(align_wd)

    # prepare scan_xyz
    coordinate1 = storage.allstruc.get_value(f"0.struc_info.general.coordinate")
    ele_l1 = storage.allstruc.get_value(f"0.struc_info.general.ele_l")
    if "1" in storage.allstruc.data.keys():
        coordinate2 = storage.allstruc.get_value(f"1.struc_info.general.coordinate")
        ele_l2 = storage.allstruc.get_value(f"1.struc_info.general.ele_l")
    else:
        coordinate2 = storage.allstruc.get_value(f"0.struc_info.general.coordinate")
        ele_l2 = storage.allstruc.get_value(f"0.struc_info.general.ele_l")

    # generate reactiv complex
    aligningcoord1, ele_l_merge = alignment(coordinate1, ele_l1, coordinate2, ele_l2, atl,
                                            extra_setting={"maxovlap": False,
                                                           "aligning_distance": storage.setting.get_value(
                                                               "searching_setting.inter.aligning_distance", 5)},
                                            writingP=True)
    write_xyz(aligningcoord1, ele_l_merge, f"coord_default")
    alignedcoord2, ele_l_merge = alignment(coordinate1, ele_l1, coordinate2, ele_l2, atl,
                                           extra_setting={"maxovlap": True,
                                                          "aligning_distance": storage.setting.get_value(
                                                              "searching_setting.inter.aligning_distance", 5)})
    write_xyz(alignedcoord2, ele_l_merge, f"coord_maxovelap")
