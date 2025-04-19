import aRST.script.initialization as storage
from aRST.script.exploration.assign_reaction import SingleCycleSearch
from aRST.script.spr import SPR
from aRST.script.reactivity.print_stuff import *
from aRST.script.toolbox import create_folder, write_xyz
from aRST.script.geom import alignment

def aligning_two_structures(reactants=None,at1=None,at2=None):
    if reactants is None:
        reactants = storage.allstruc.reactants_ini
    if len(reactants) !=2:
        raise ValueError("The number of input structures is not 2, please check")
    if at1 is None and at2 is None:
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

        react_l = react_l[0]
        # react_l is ({"struc_combi": [str(istruc_l), str(istruc_h)],
        #                 "at_combi": [[str(iat_l), str(iat_h)]],
        #                 "HOAO": H,
        #                 "LUAO": L,
        #                 "gap": gap})
        if int(react_l["struc_combi"][0]) == 0 :
            at1 = int(react_l["at_combi"][0][0])
            at2 = int(react_l["at_combi"][0][1])
        else:
            at1 = int(react_l["at_combi"][0][1])
            at2 = int(react_l["at_combi"][0][0])
    else:
        print(f"pre-assigned reactive atom pair is: struc 0 : {at1}, struc 1: {at2}")

    align_wd = join(storage.setting.wd0, f"alignment")
    create_folder(align_wd)

    # prepare scan_xyz
    coord1 = storage.allstruc.get_value(f"0.struc_info.general.coord")
    coord2 = storage.allstruc.get_value(f"1.struc_info.general.coord")

    # generate reactiv complex
    aligningcoord1 = alignment(coord1, coord2, at1, at2,
                             extra_setting={"maxovlap": False,
                                 "aligning_distance": storage.setting.get_value(
                                     "searching_setting.inter.aligning_distance", 5)})
    write_xyz(aligningcoord1, f"coord_default")
    alignedcoord2 = alignment(coord1, coord2, at1, at2,
                             extra_setting={"maxovlap": True,
                                 "aligning_distance": storage.setting.get_value(
                                     "searching_setting.inter.aligning_distance", 5)})
    write_xyz(alignedcoord2, f"coord_maxovelap")
