import toml
import threading
import os
from os.path import join, exists
from aRST.script.toolbox import  read_coord
from aRST.script.callsys import  CallMolBar
from aRST.script.exploration.structure import Allstruc,Singlestruc
from aRST.script import input_class
from pydantic import ValidationError


class Setting:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, infile=None):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super().__new__(cls)
                cls._instance._load_config(infile)
        return cls._instance

    def _load_config(self, infile):
        with open(infile, "r") as f:
            self.config = toml.load(f)

        # checking config
        print("Checking input TOML...", end="")
        try:
            self.model = input_class.InputModel.model_validate(self.config)
        except ValidationError as e:
            print("Input file validation failed:")
            print(e)
            raise SystemExit(1)
        print("done.\n")

        # checking env
        print("Checking envirment:")
        self.check_environment()
        print("done.\n")

        # checking the validty of old restart files
        check_restart_compatible(self.config, self.model.searching_setting.general.wd)

        # clean up buffer struc
        if exists(self.model.searching_setting.general.bufferwd):
            import shutil
            shutil.rmtree(self.model.searching_setting.general.bufferwd)

    def check_environment(self):
        # xtb is necessary
        xtb_bin = os.environ.get("XTB_BIN")
        if not xtb_bin:
            raise RuntimeError("XTB_BIN is not set.")
        print(f"XTB_BIN checked, PATH={xtb_bin}")

        # stda is necessary
        stda_bin = os.environ.get("STDA_BIN")
        if not stda_bin:
            raise RuntimeError("STDA_BIN is not set.")
        print(f"stda_bin checked, PATH={stda_bin}")

        # gxtb
        if self.model.searching_setting.general.sp_method == "gxtb":
            gxtb_bin = os.environ.get("GXTB_BIN")
            if not gxtb_bin:
                raise RuntimeError("sp_method is 'gxtb' but GXTB_BIN is not set.")
            print(f"GXTB_BIN checked, PATH={gxtb_bin}")

        reactivity = self.model.searching_setting.reactivity
        # orca
        if reactivity.MO_method == input_class.MOMethod.orca:
            orca_bin = os.environ.get("ORCA_BIN")
            if not orca_bin:
                raise RuntimeError(
                    "MO_method is 'orca' but environment variable ORCA_BIN is not set."
                )
            print(f"ORCA_BIN checked, PATH={orca_bin}")
            mwfn_bin = os.environ.get("MWFN_BIN")
            if not mwfn_bin:
                raise RuntimeError(
                    "MO_method with 'orca' needs Multiwfn for converting basis functions but MWFN_BIN is not set."
                )
            print(f"MWFN_BIN checked (!! make sure you have set iloadasCart= 1 !!), PATH={mwfn_bin}")
        # terachem
        elif reactivity.MO_method == input_class.MOMethod.terachem:
            tc_bin = os.environ.get("TC_BIN")
            if not tc_bin:
                raise RuntimeError(
                    "MO_method is 'terachem' but environment variable TC_BIN is not set."
                )
            print(f"TC_BIN checked, PATH={tc_bin}")

        # crest
        disso = self.model.searching_setting.disso
        if disso and (disso.depotonation or disso.potonation
        ):
            crest_bin = os.environ.get("CREST_BIN")
            if not crest_bin:
                raise RuntimeError(
                    "(de)potonation set but environment variable CREST_BIN is not set."
                )
            print(f"CREST_BIN checked, PATH={crest_bin}")



def register_init_input_reactants(setting,allstruc,getfeatures=True):
    reactants_ini = []  # record reactants for the input
    wd = setting.searching_setting.general.wd
    ## reactants
    reactant_keys = [k for k in setting.geom.keys() if k.startswith("reactant")]
    for id, key in enumerate(reactant_keys):
        reactant = setting.geom[key]
        coord_path = join(wd, reactant.name)
        coordinate, ele_l = read_coord(coord_path)
        c = reactant.charge
        m = reactant.multiplicity
        mb = CallMolBar(coordinate, ele_l).get_mb()
        ns = setting.geom[key].number  # repeat adding time

        currentstruc = Singlestruc(coordinate,ele_l,c,m,mb)

        sid_exist, refsid = currentstruc.structure_exist()
        if not sid_exist:
            sid,s = allstruc.placeholde_new_sid()
            if setting.read.general.jumprefine:
                features,readsuss = currentstruc.read_features_from_struc_folder(rootfolder=join(wd, "struc", str(sid)),
                                                                                 charge=c,m=m)
                if readsuss:
                    allstruc.register_new_sid(sid,s, features)
                else:
                    if getfeatures:
                        features = currentstruc.process_workflow(rootfolder=join(wd,"struc",str(sid)),ideal=False)
                    else:
                        features = {
                            "coordinate": coordinate,
                            "ele_l": ele_l,
                            "multiplicity": m,
                            "charge": c,
                            "mb": mb,
                        }
            else:
                if getfeatures:
                    features = currentstruc.process_workflow(rootfolder=join(wd, "struc", str(sid)),ideal=False)
                else:
                    features = {
                        "coordinate": coordinate,
                        "ele_l": ele_l,
                        "multiplicity": m,
                        "charge": c,
                        "mb": mb,
                    }
            allstruc.register_new_sid(sid,s, features)
            for _ in range(ns):
                reactants_ini.append(sid)
        else:
            for _ in range(ns):
                reactants_ini.append(refsid)
    return reactants_ini


setting = None
allstruc = None

import hashlib
import json
import shutil
def make_config_hash(config):
    EXCLUDE_KEYS = {
        ("searching_setting", "disso", "number"),
        ("searching_setting", "inter", "number"),
        ("searching_setting", "intra", "number"),
    }
    filtered = filter_config(config, EXCLUDE_KEYS)

    import json, hashlib
    text = json.dumps(filtered, sort_keys=True, default=str)
    return hashlib.sha256(text.encode()).hexdigest()


def filter_config(config, exclude_keys):
    def _filter(d, path=()):
        if not isinstance(d, dict):
            return d

        new_d = {}
        for k, v in d.items():
            new_path = path + (k,)

            if new_path in exclude_keys:
                continue

            new_d[k] = _filter(v, new_path)

        return new_d

    return _filter(config)

def check_restart_compatible(config: dict, wd: str, restart_dir="checkpoint"):
    restart_path = join(wd, restart_dir)
    hash_file = join(restart_path, "input.sha256")

    new_hash = make_config_hash(config)

    if exists(hash_file):
        with open(hash_file, "r") as f:
            old_hash = f.read().strip()

        if old_hash != new_hash:
            print("NOTICE: input TOML changed, removing old restart files")
            shutil.rmtree(restart_path)

    os.makedirs(restart_path, exist_ok=True)

    with open(hash_file, "w") as f:
        f.write(new_hash)

def initialize(setting_instance):
    """Initialize global setting and allstruc with the provided Setting instance."""
    global setting, allstruc
    setting = setting_instance.model  # Assign the passed Setting instance
    allstruc = Allstruc(setting)