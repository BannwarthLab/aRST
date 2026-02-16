from aRST.script.exploration.assign_reaction import SingleCycleSearch
from aRST.script.reactivity.print_stuff import *
from aRST.script.exploration.interA import *
from aRST.script.exploration.intraA import *
from aRST.script.exploration.intraD import *
def show_afo(printnumber=10):
    reactants = storage.allstruc.reactants_ini
    # preapre afo
    SingleCycleSearch(reactants)
    coordinates = [storage.allstruc.get_value(f"{_}.struc_info.general.coordinate") for _ in reactants]
    ele_ls = [storage.allstruc.get_value(f"{_}.struc_info.general.ele_l") for _ in reactants]

    # show interA reactivity
    if len(reactants) > 1:
        react_l = SPR(reactants, storage.setting.get_value("searching_setting.general.spr_mode", 0)).inter(
            equivalent_afogap_tol=storage.setting.get_value("searching_setting.general.equivalent_afogap_tol", 0.01),
            exclude_elements=storage.setting.get_value("searching_setting.general.exclude_elements"),
            exclude_pairs=storage.setting.get_value("searching_setting.general.exclude_pairs")
        )
        # write full react_l in file
        with open(join(storage.setting.wd0, "react_l_interA"), "w") as f:
            for item in react_l:
                f.write(json.dumps(item) + "\n")

    # assess intraA reactivity
    react_l = SPR(reactants,
                  storage.setting.get_value("searching_setting.general.spr_mode", 0)
                  ).intra(
        equivalent_afogap_tol=storage.setting.get_value("searching_setting.general.equivalent_afogap_tol", 0.01),
        exclude_elements=storage.setting.get_value("searching_setting.general.exclude_elements"),
        exclude_pairs=storage.setting.get_value("searching_setting.general.exclude_pairs")
    )
    # write full react_l in file
    with open(join(storage.setting.wd0, "react_l_intraA"), "w") as f:
        for item in react_l:
            f.write(json.dumps(item) + "\n")

    # assess Disso reactivity
    react_l = SDP(reactants,
                  exclude_elements=storage.setting.get_value("searching_setting.general.exclude_elements"),
                  exclude_pairs=storage.setting.get_value("searching_setting.general.exclude_pairs")
                  ).build_react_l(tol=0.01)
    # write full react_l in file
    with open(join(storage.setting.wd0, "react_l_disso"), "w") as f:
        for item in react_l:
            f.write(json.dumps(item) + "\n")
