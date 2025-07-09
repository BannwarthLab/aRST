import sys
from aRST.script.toolbox import mydict, get_all_keys
import os
import shutil

required_input = {
    "read": {
        "general": {
            "dftopt": (bool, False),  # doing dftopt when first reading input coord
            "gfn2opt": (bool, False),  # doing gfn2 when first reading input coord
            "jumprefine": (bool, False),  # direct read struc info from folder struc/0,1,2,.....
        }
    },
    "searching_setting": {
        "general": {
            "cycle": (int, 1),  # reaction simulation cycle
            "elim": ((int, float), 300),
            # energy limitation to kill this reaction path in current path, unit in kcal/mol
            "hydrogen": (bool, True),  # include hydrogen element in the AFO searching
            "halogen": (bool, True),  # include halogen elements in the AFO searching
            "spr_mode": (int, 0),
            # AFO gap formation rules: 0: take both spin same/oppo combinations into account; 1: only same spin; 2: only oppo spin
            "scan_mode": (int, 0),  # bond formation method: 0: using xtb scan; 1: using harmonic FF
            "sep_mode": (int, 0),
            # separation rules for preliminary product structure: 0: sep fragments before gfn2 opt; 1: opt product before separation
            "mb_level": (str, None),  # "loose" to only screen with topology spectrum
            "keepmolden": (bool, False)  # del all molden files to avoid large I/O writing
        },
        "reactivity": {
            "afo_method": (str, "gfn2"),  # method for afo feature calculations
        },
        "inter": {
            "number": (int, 0),  # number of intermolecular reaction paths to be searched in one cycle
            "search_mode": (int, 1),
            # searching rules: 0: only perform bonding simulation to same structure (A+A); 1: ..to different structure (A+B); 2: both (A+A and A+B)
            "at1": (int, None),
            # manually define the search pair, index strarts from 0. at1,at3 belong to struc 0; at2,at4 belong to struc1
            "at2": (int, None),
            "at3": (int, None),
            "at4": (int, None),
            "loose_bond": (bool, False),  # loose a bit the weakest bond in the structure after traget bond formation
            "maxovlap": (bool, False),  # alignment method for achiving maximum overlap of reacting structures
            "aligning_distance": ((int, float), 5),  # scale of aligning distance of reacting structures
            "concerted_pair": (bool, False)  # reate concerted pair in simulations (good for DA reactions)
        },
        "intra": {
            "number": (int, 0),  # number of intramolecular reaction paths to be searched in one cycle
            "at1": (int, None),
            # manually define the search pair, index strarts from 0. at1,at3 belong to struc 0; at2,at4 belong to struc1
            "at2": (int, None),
            "at3": (int, None),
            "at4": (int, None),
            "concerted_pair": (bool, False)  # reate concerted pair in simulations (good for DA reactions)
        },
        "disso": {
            "number": (int, 0),  # number of dissociative reaction paths to be searched in one cycle
            "at1": (int, None),
            # manually define the search pair, index strarts from 0. at1,at3 belong to struc 0; at2,at4 belong to struc1
            "at2": (int, None),
            "at3": (int, None),
            "at4": (int, None),
            "mode": (int, 0),  # 0: use mbo as bo; 1: use embo as bo for ranking bond term strength
            "concerted_pair": (bool, False),  # reate concerted pair in simulations (good for DA reactions)
            "pulling_distance": ((int, float, str), "auto"),  # scale of pulling distance
            "depotonation": (bool, False),  # checking depotonated reactant before dissociation
            "potonation": (bool, False),  # checking potonated reactant before dissociation
        }
    },
    "xtb_setting": {  # this settings go into both gfn sp and opt calculations
        "method": (str, "gfn2"),  # GFN lvl
        "optlevel": (str, None),  # xtb acceptable optimization levels
        "alpb": (str, "ch2cl2"),  # ALPB solvation model, aRST only pass keyword to xtb program
        "struc_etemp": (int, 300),  # electronic temp for structure opt, recommended using RT temp.
    },
    "orca_setting": {  # this settings go into both DFT sp and opt calculations
        "functional": (str, "r2scan-3c"),  # functional method
        "basis": (str, None),  # basis set
        "dftopt": (bool, False),  # if doing DFT opt for product structure
        "solmode": (str, "ch2cl2"),  # CPCM solvation model, aRST only pass keyword to orca program
        "ncpu": (int, 1),  # parallel calculations for orca
    }
}

condition_md = [
    {
        "condition": [("searching_setting", "general", "MD", True)],
        "required_block": {
            "searching_setting": {
                "MD": (bool, False),  # perform short md simulations to test the stability of product structure
                "md": {
                    "time": (int, 10),  # total run time of MD simulation, unit in ps
                    "step": (int, 1),  # time step for propagation
                    "temp": (int, 300),  # electronic temp, by defalut using scan temp, unit in K
                    "shake": (int, 1),  # 1: constrain H bonds; 2: constrain all bonds
                    "dump": (int, 50),  # interval for trajectory printout
                }
            }
        },
        "message": "Since MD simulation is asked, the block [searching_setting.md] must be provided."
    }
]
condition_AFO_orca = [
    {
        "condition": [("searching_setting", "reactivity", "afo_method", "stda"),
                      ("searching_setting", "reactivity", "MO_method", "orca")],
        "required_block": {
            "searching_setting": {
                "reactivity": {
                    "afo_method": (str, "stda"),
                    "MO_method": (str, "orca"),
                    "afo_orca_functional": (str, "bhlyp"),  # functional method
                    "afo_orca_basis": (str, None),  # basis set
                    "afo_orca_cpcm": (str, "ch2cl2"),  # CPCM solvation model, aRST only pass keyword to orca program
                }
            }
        },
        "message": "Since AFO features are asked to calculated at DFT methods in ORCA, the DFT methods must be provided."
    }
]
condition_AFO_terachem = [
    {
        "condition": [("searching_setting", "reactivity", "afo_method", "stda"),
                      ("searching_setting", "reactivity", "MO_method", "terachem")],
        "required_block": {
            "searching_setting": {
                "reactivity": {
                    "afo_method": (str, "stda"),
                    "MO_method": (str, "terachem"),
                    "afo_orca_functional": (str, "bhlyp"),  # functional method
                    "afo_orca_basis": (str, None),  # basis set
                    "afo_tc_epsilon": ((float, int), 9.08),
                    # CPCM solvation model, aRST only pass keyword to Terachem program
                }
            }
        },
        "message": "Since AFO features are asked to calculated at DFT methods in Terachem, the DFT methods must be provided."
    }
]
condition_xtbscan = [
    {
        "condition": [("searching_setting", "general", "scan_mode", 0)],
        "required_block": {
            "xtb_setting": {
                "scan_etemp": (int, 1000),  # electronic temp for xtbscan/harmonic FF, recommended using higher temp.
                "k": (int, 1),  # geom constrain strength
                "step": (int, 20),  # pushing/pulling steps
            }
        },
        "message": ""
    }
]
condition_FF = [
    {
        "condition": [("searching_setting", "general", "scan_mode", 1)],
        "required_block": {
            "xtb_setting": {
                "scan_etemp": (int, 1000),  # electronic temp for xtbscan/harmonic FF, recommended using higher temp.
                "k": (int, 1),  # geom constrain strength
            }
        },
        "message": ""
    }
]

conditions = []
for i in (condition_FF, condition_xtbscan, condition_md, condition_AFO_terachem, condition_AFO_orca):
    conditions.append(i)


def checking_parameter(setting):
    assert isinstance(setting, dict) and setting
    schema = mydict(required_input)
    # config = mydict(setting)
    for item in conditions:
        item = item[0]
        condition = item["condition"]
        add = True
        for item2 in condition:
            key = ".".join(list(item2[:-1]))
            expect_val = item2[-1]
            if isinstance(expect_val, str):
                expect_val = expect_val.lower()
            if setting.get_value(key) != expect_val:
                add = False
                break
        if add:
            for block_key, block_value in item["required_block"].items():
                if block_key not in required_input:
                    required_input[block_key] = {}
                for subkey, subval in block_value.items():
                    required_input[block_key][subkey] = subval

    all_keys_in_config = get_all_keys(setting)
    all_keys_in_shema = get_all_keys(schema)
    allcorrect = True
    for item in all_keys_in_config:
        # check if input key exists
        if item in all_keys_in_shema:
            # check if type is correct
            val = setting.get_value(item)
            expected_type = schema.get_value(item)[0]
            if isinstance(expected_type, (list, tuple)):
                accepted_type = list(expected_type)
            else:
                accepted_type = [expected_type]

            if type(val) in accepted_type:
                pass
            else:
                allcorrect = False
                print("Error type in", val, ": This type is not accepted for", item, f"{accepted_type}")
        elif "geom" in item:
            # geom block will be checked in init
            continue
        elif "wd" in item:
            # accepted in aRST 1.0
            continue
        else:
            allcorrect = False
            print("Error parameter in ", item, ": This parameter is not accept in aRST")
    if not allcorrect:
        raise ValueError("wrong input!")
    else:
        # adding missing parameter to default val
        missing_par = [x for x in all_keys_in_shema if x not in all_keys_in_config]
        print("The following missing parameters will be set to their default values:")
        for par in missing_par:
            defalut_val = schema.get_value(par)[1]
            if defalut_val is not None and not (isinstance(defalut_val, bool) and not defalut_val):
                setting.set_value(par, defalut_val)
                print(f"{par} -> {defalut_val}")
        print()
    return setting



def checking_env(setting):
    required_vars = ["ORCA_BIN","XTB_BIN"]#,"TRVT_BIN"
    if setting.get_value("searching_setting.reactivity.afo_method") =="stda":
        required_vars.append("STDA_BIN")
        required_vars.append("MWFN_BIN")
        if setting.get_value("searching_setting.reactivity.MO_method") =="terachem":
            required_vars.append("TC_BIN")
    elif setting.get_value("searching_setting.reactivity.afo_method") =="gfn2":
        required_vars.append("XTBML_BIN")
    if setting.get_value("searching_setting.disso.depotonation") or setting.get_value("searching_setting.disso.potonation"):
        required_vars.append("CREST_BIN")

    missing = [var for var in required_vars if var not in os.environ]
    if missing:
        print("Missing environment variables:", ", ".join(missing))
        print("Please set them using 'export VAR=VALUE' in your shell.")
        exit(1)

