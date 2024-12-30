import os.path
from aRST.script.toolbox import create_folder
from aRST.script.dictionaryInspector import DictionaryInspector
from aRST.script.assign_reaction import SingleCycleSearch
from aRST.script.spr import SPR_HL, SPR_dissociation,build_bo_d,find_identical_at
import json


def get_reactivity(wd, STRUC_rec, setting_d, reactants):
    # borrow ml prepare from singlecycle class
    singlecycle = SingleCycleSearch(STRUC_rec, {}, setting_d, wd)
    STRUC_inspector = DictionaryInspector(STRUC_rec)
    setting_inspector = DictionaryInspector(setting_d)
    for struc in reactants:
        if STRUC_inspector.get_value(f"{struc}.struc_info.reactivity.ml") is None or STRUC_inspector.get_value(
                f"{struc}.struc_info.reactivity.bo") is None:
            singlecycle.ML_prepare(struc)

    if len(reactants) > 1:
        # get interA reaction possibility
        react_l = SPR_HL(reactants,
                         STRUC_rec,
                         allow_hydrogen=setting_inspector.get_value("searching_setting.inter.hydrogen", True),
                         allow_halogen=setting_inspector.get_value("searching_setting.inter.halogen", True),
                         search_mode=setting_inspector.get_value("searching_setting.inter.search_mode", 1),
                         spr_mode=setting_inspector.get_value("searching_setting.inter.spr_mode", 0)
                         )



        print("starting to search inter-association possibility:")
        print("allow_hydrogen:",setting_inspector.get_value("searching_setting.inter.hydrogen", True))
        print("allow_halogen:", setting_inspector.get_value("searching_setting.inter.halogen", True))
        print("search_range:", setting_inspector.get_value("searching_setting.inter.number", 1))
        print("write jason and txt file in ", os.path.join("reactivity", "interA"))
        print("\n")
        if react_l:
            with open(os.path.join(wd, "reactivity_interA"), "w") as f:
                for entry in react_l:
                    f.write(json.dumps(entry))
                    f.write('\n')
        else:
            raise Exception(f"no possible interA possibilities can be found in reactants list: {reactants}.")

    if setting_inspector.get_value("searching_setting.disso"):
        print("starting to search dissociation possibility:")
        print("allow_hydrogen:", setting_inspector.get_value("searching_setting.disso.hydrogen", True))
        print("allow_halogen:", setting_inspector.get_value("searching_setting.disso.halogen", True))
        print("search_range:", setting_inspector.get_value("searching_setting.disso.number", 1))
        print("search_mode:", setting_inspector.get_value("searching_setting.disso.search_mode", "bo"))
        print("write jason and txt file in ", os.path.join("reactivity", "disso"))
        print("for jason dict: {BO:[strucid, at1, at2 , dgap]}")
        print("\n")
        BO = SPR_dissociation(reactants,
                              STRUC_rec,
                              allow_hydrogen=setting_inspector.get_value(
                                  "searching_setting.disso.hydrogen", True),
                              allow_halogen=setting_inspector.get_value(
                                  "searching_setting.disso.halogen", True),
                              search_mode=0 if str.lower(
                                  setting_inspector.get_value("searching_setting.disso.mode", "wboe")) == "wboe" else 1)
        if BO:
            outputwd = create_folder(os.path.join(wd, "reactivity", "disso"))
            with open(f'{outputwd}/disso.json', "w") as f:
                f.write(json.dumps(BO, indent=2))
                f.write('\n')

            #[strucid, at1 - 1, at2 - 1, dgap]
            with open(f'{outputwd}/top10', "w") as f:
                f.write("    gap,    STRUC,    AT1,    AT2, degenerated info\n")
                for key in sorted(BO.keys())[:setting_inspector.get_value("searching_setting.disso.number", 1)]:
                    if len(BO[key])==1:
                        struc,at1,at2 = BO[key][0][:3]
                        f.write(f"%8.3f,%6d,%6d,%6d\n" % (key,int(struc),int(at1), int(at2)))
                    else:
                        f.write(f"%8.3f,"%key)
                        for value in BO[key]:
                            struc, at1, at2 = value[:3]
                            f.write(f"%6d,%6d,%6d; " % (int(struc), int(at1), int(at2)))
                        f.write("\n")
        else:
            raise Exception(f"no possible dissociation possibilities can be found in reactants list: {reactants}.")

    if setting_inspector.get_value("searching_setting.intra"):
        react_aa, react_bb, react_ab, react_ba = SPR_HL(reactants,
                                                        STRUC_rec,
                                                        allow_hydrogen=setting_inspector.get_value(
                                                            "searching_setting.intra.hydrogen", True),
                                                        allow_halogen=setting_inspector.get_value(
                                                            "searching_setting.intra.halogen", True),
                                                        search_range=0
                                                        )
        bo = build_bo_d(reactants,
                        STRUC_inspector,
                        allow_hydrogen=setting_inspector.get_value(
                            "searching_setting.intra.hydrogen", True),
                        allow_halogen=setting_inspector.get_value(
                            "searching_setting.intra.halogen", True))
        print("starting to search intra association possibility:")
        print("allow_hydrogen:", setting_inspector.get_value("searching_setting.intra.hydrogen", True))
        print("allow_halogen:", setting_inspector.get_value("searching_setting.intra.halogen", True))
        print("search_range:", setting_inspector.get_value("searching_setting.intra.number", 1))
        print("write jason and txt file in ", os.path.join("reactivity", "intraA"))
        print("for jason dict: "
              "{'hoao-luao gap':[struc index1, struc index2,[[at index 1,at index2],[possible eq atomic combination]]]")
        print("\n")
        if any([react_aa, react_bb, react_ab, react_ba]) and bo:
            outputwd = create_folder(os.path.join(wd, "reactivity", "intraA"))
            name_l = ["aa", "bb", "ab", "ba"]
            for react_d, name in zip([react_aa, react_bb, react_ab, react_ba], name_l):
                if react_d:
                    intraa = {}
                    key_l = sorted(react_d.keys())[:setting_inspector.get_value(
                        "searching_setting.intra.number", 10)]
                    for current_key in key_l:
                        current_item = react_d[current_key]
                        for struc_combi, at_l in current_item.items():
                            struc = struc_combi[0]
                            if 1 not in [len(set(k)) for k in at_l]:  # not same atom
                                for pair in at_l:
                                    if pair in bo[struc] or pair[::-1] in bo[struc]:  # not bonded
                                        continue
                                    else:
                                        intraa[current_key] = struc, at_l
                    with open(f'{outputwd}/{name}.json', "w") as f:
                        f.write(json.dumps(intraa, indent=2))
                        f.write('\n')
                    with open(f'{outputwd}/top_{name}', "w") as f:
                        f.write("    gap,    STRUC1,     at_l\n")
                        for gap in sorted(intraa.keys())[:10]:
                            struc = intraa[gap][0]
                            at_l = intraa[gap][1]
                            f.write(f"%8.3f,%6d,      {at_l}\n" % (gap, int(struc)))

    return







