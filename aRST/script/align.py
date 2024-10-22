import os
from aRST.script.dictionaryInspector import DictionaryInspector
from aRST.script.toolbox import create_folder, write_xyz
from aRST.script.geom import alignment


def aligning(wd, STRUC_rec, setting_d, reactants):
    outputwd = create_folder(os.path.join(wd, "align"))
    struc_inspector = DictionaryInspector(STRUC_rec)
    setting_inspector = DictionaryInspector(setting_d)

    struc1, struc2 = reactants
    at1 = struc_inspector.get_value(f"{struc1}.struc_info.general.at")
    at2 = struc_inspector.get_value(f"{struc2}.struc_info.general.at")

    coord1 = struc_inspector.get_value(f"{struc1}.struc_info.general.coord")
    coord2 = struc_inspector.get_value(f"{struc2}.struc_info.general.coord")

    alignedcoord = alignment(coord1, coord2, at1, at2,
                             extra_setting={"maxovlap": setting_inspector.get_value(
                                 "searching_setting.inter.maxovlap", False),
                                 "aligning_distance": setting_inspector.get_value(
                                     "searching_setting.inter.aligning_distance", 5)},writingP=True)
    write_xyz(alignedcoord, f"aligned_coord")
