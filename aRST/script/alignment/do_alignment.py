import aRST.script.initialization as storage
from aRST.script.exploration.reactionsearch import SingleCycleSearch
from aRST.script.spr import SPR
from aRST.script.reactivity.print_stuff import *
from aRST.script.toolbox import create_folder, write_xyz
from aRST.script.geom import alignment


def aligning_two_structures(reactants, atl=None):
    from aRST.script.exploration.structure import structure_class
    if len(reactants) != 2:
        raise ValueError("The number of input structures is not 2, please check")
    if atl is None:
        # check if define reactive ats
        atom_pairs = storage.setting.searching_setting.inter.atom_pair
        if atom_pairs:
            # overwrite react_l with user's input
            # predifined reaction, used in one cycle
            at1, at2 = atom_pairs[0]
            react_l = [{"reactants": [reactants[0], reactants[1]],
                        "atl": [[at1, at2]],
                        "gap": 0,
                        "force_combination": True}]
            if len(atom_pairs) > 1:
                storage.setting.searching_setting.inter.concerted_pair=True
                for atom_pair in atom_pairs[1:]:
                    at_combi_list = react_l[0]["atl"]
                    at3, at4 = atom_pair
                    at_combi_list.append([at3, at4])
        else:
            print("aRST doesn't detect preassign reactive ats, switch to most reactive pair assigning by AFO gaps")
            SingleCycleSearch(reactants)
            spr = SPR(reactants, storage.setting.searching_setting.general.spr_mode)
            react_l = spr.inter(
                equivalent_afogap_tol=storage.setting.searching_setting.general.equivalent_afogap_tol,
                exclude_elements=storage.setting.searching_setting.general.exclude_elements,
                exclude_pairs=storage.setting.searching_setting.general.exclude_pairs)
            cut = 2 if storage.setting.searching_setting.inter.concerted_pair else 1
            for item in react_l:
                pairs = item["atl"][:cut]
                item["atl"] = pairs
            entry = react_l[0]
            at_combi_list = entry["atl"]

    align_wd = join(storage.setting.searching_setting.general.wd, f"alignment")
    create_folder(align_wd)

    # prepare scan_xyz
    s1 = storage.allstruc.get(0)
    coordinate1 = s1.general.coordinate
    ele_l1 = s1.general.ele_l

    if 1 in storage.allstruc.data.keys():
        s2 = storage.allstruc.get(1)
    else:
        s2 = storage.allstruc.get(0)

    coordinate2 = s2.general.coordinate
    ele_l2 = s2.general.ele_l


    # generate reactiv complex
    aligningcoord1, ele_l_merge = alignment(coordinate1, ele_l1, coordinate2, ele_l2, at_combi_list,
                                            extra_setting={"maxovlap": False,
                                                           "aligning_distance": storage.setting.searching_setting.inter.aligning_distance},
                                            writingP=True)
    write_xyz(aligningcoord1, ele_l_merge, f"coord_default")
    alignedcoord2, ele_l_merge = alignment(coordinate1, ele_l1, coordinate2, ele_l2, at_combi_list,
                                           extra_setting={"maxovlap": True,
                                                          "aligning_distance": storage.setting.searching_setting.inter.aligning_distance})
    write_xyz(alignedcoord2, ele_l_merge, f"coord_maxovelap")
