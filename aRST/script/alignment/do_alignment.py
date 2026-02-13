import aRST.script.initialization as storage
from aRST.script.exploration.assign_reaction import SingleCycleSearch
from aRST.script.spr import SPR
from aRST.script.reactivity.print_stuff import *
from aRST.script.toolbox import create_folder, write_xyz
from aRST.script.geom import alignment

def aligning_two_structures(reactants=None,atl=None):
    if reactants is None:
        reactants = storage.allstruc.reactants_ini
    if len(reactants) !=2:
        raise ValueError("The number of input structures is not 2, please check")
    if atl is None:
        # check if define reactive ats
        at1 = storage.setting.get_value("searching_setting.inter.at1")
        at2 = storage.setting.get_value("searching_setting.inter.at2")
        if at1 is None or at2 is None:
            print("aRST doesn't detect preassign reactive ats, switch to most reactive pair assigning by AFO gaps")
            SingleCycleSearch(reactants)
            hydrogen = storage.setting.get_value("searching_setting.general.hydrogen", True)
            halogen = storage.setting.get_value("searching_setting.general.halogen", True)
            print(f"allow hydrogen atoms: {hydrogen}\nallow halogen atoms: {halogen}")
            spr = SPR(reactants, hydrogen, halogen)
            react_l = spr.inter(storage.setting.get_value("searching_setting.general.spr_mode", 0))
            cut = 2 if storage.setting.get_value("searching_setting.inter.concerted_pair", False) else 1
            for item in react_l:
                pairs = item["at_combi"][:cut]
                item["at_combi"] = pairs
            entry = react_l[0]
            atl = entry["at_combi"]
        else:
            atl = [[int(at1),int(at2)]]
            at3 = storage.setting.get_value("searching_setting.inter.at3")
            at4 = storage.setting.get_value("searching_setting.inter.at4")
            if at3 is not None and at4 is not None:
                atl.append([int(at3),int(at4)])

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
    aligningcoord1,ele_l_merge = alignment(coordinate1,ele_l1,coordinate2,ele_l2, atl,
                             extra_setting={"maxovlap": False,
                                 "aligning_distance": storage.setting.get_value(
                                     "searching_setting.inter.aligning_distance", 5)},writingP=True)
    write_xyz(aligningcoord1,ele_l_merge, f"coord_default")
    alignedcoord2,ele_l_merge = alignment(coordinate1,ele_l1,coordinate2,ele_l2, atl,
                             extra_setting={"maxovlap": True,
                                 "aligning_distance": storage.setting.get_value(
                                     "searching_setting.inter.aligning_distance", 5)})
    write_xyz(alignedcoord2,ele_l_merge, f"coord_maxovelap")
