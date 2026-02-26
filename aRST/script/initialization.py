import toml
import threading
import os
from typing import Dict
from os.path import join, exists
from aRST.script.toolbox import create_folder, mydict, read_coord, get_strucinfo_from_path, write_xyz
from aRST.script.callsys import CallxTB, CallORCA, CallgXTB, CallMolBar
from aRST.script import input_class, structure_class
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
        self.config = mydict(self.config)

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


class Allstruc:
    def __init__(self, model):
        self.model = model
        self.data = mydict({})
        self.reactants_ini = []  # record reactants for the input
        self._next_id = 0
        self._lock = threading.Lock()
        self.mb_index: Dict[str, int] = {}  # <-- dedup index
        self._initialize()


    def create_new_structure(self):
        with self._lock:
            sid = self._next_id
            self._next_id += 1

        s = structure_class.Structure(id=sid)
        self.data[sid] = s
        return sid,s


    def get_or_create_new_structure(self, mb):
        # return data, id, created
        with self._lock:
            sid = self.mb_index.get(mb)
            if sid is not None:
                return self.data[sid], sid, False  # existing

            sid = self._next_id
            self._next_id += 1
            s = structure_class.Structure(id=sid)
            self.data[sid] = s
            self.mb_index[mb] = sid
            return s, sid, True  # newly created

    def _initialize(self):
        wd = self.model.searching_setting.general.wd
        ## reactants
        reactant_keys = [k for k in self.model.geom.keys() if k.startswith("reactant")]
        for id, key in enumerate(reactant_keys):
            reactant = self.model.geom[key]
            coord_path = join(wd, reactant.name)
            coordinate, ele_l = read_coord(coord_path)
            c = reactant.charge
            m = reactant.multiplicity
            mb = CallMolBar(coordinate, ele_l).get_mb()

            # set up new id if needed
            s, sid, created = self.get_or_create_new_structure(mb)
            if created:
                s.struc_info.general.coordinate = coordinate
                s.struc_info.general.ele_l = ele_l
                s.struc_info.general.charge = c
                s.struc_info.general.multiplicity = m
                s.struc_info.general.mb = mb
            self.reactants_ini.append(sid)
            ns = self.model.geom[key].number
            s.struc_info.general.number = ns
            for _ in range(ns-1):
                self.reactants_ini.append(sid)

        if self.model.read.general.jumprefine:
            # collect struc refine info
            print(f"Reading structure information of input {list(set(self.reactants_ini))} from folder struc/ ...")
            for sid in list(set(self.reactants_ini)):
                struc_info, info_complete = get_strucinfo_from_path(join(wd, "struc", str(sid)))
                if info_complete:
                    # compare mb to make sure
                    if struc_info["mb"] != self.data[sid].struc_info.general.mb:
                        print(self.data[sid].struc_info.general.mb)
                    self.set_value_from_struc_info(sid, struc_info)
                else:
                    print(f"Info of structure {sid} is not complete, start to refine it")
                    coordinate = self.data[sid].struc_info.general.coordinate
                    ele_l = self.data[sid].struc_info.general.ele_l
                    c = self.data[sid].struc_info.general.charge
                    m = self.data[sid].struc_info.general.multiplicity
                    struc_info = self.refine_struc_info(sid, coordinate, ele_l, c, m)
                    self.set_value_from_struc_info(sid, struc_info)
            print("done.\n")
        else:
            print(f"Refining structure information of input {list(set(self.reactants_ini))} ...")
            for sid in list(set(self.reactants_ini)):
                coordinate = self.data[sid].struc_info.general.coordinate
                ele_l = self.data[sid].struc_info.general.ele_l
                c = self.data[sid].struc_info.general.charge
                m = self.data[sid].struc_info.general.multiplicity
                struc_info = self.refine_struc_info(sid, coordinate, ele_l, c, m)
                self.set_value_from_struc_info(sid, struc_info)
            print("done.\n")

    def set_value_from_struc_info(self, sid, struc_info):
        self.data[sid].struc_info.general.coordinate = struc_info["coordinate"]
        self.data[sid].struc_info.general.ele_l = struc_info["ele_l"]
        self.data[sid].struc_info.general.nel = struc_info["nel"]
        self.data[sid].struc_info.general.multiplicity = struc_info["m"]
        self.data[sid].struc_info.general.charge = struc_info["c"]
        self.data[sid].struc_info.general.atcharge = struc_info["atcharge"]
        self.data[sid].struc_info.general.mb = struc_info["mb"]
        self.data[sid].struc_info.general.atomic_connection = struc_info["constrain_l"]
        self.data[sid].struc_info.energy.sp = struc_info["Esp"]

    def refine_struc_info(self, sid, coordinate, ele_l, charge, m):
        comdir = join(self.model.searching_setting.general.wd, "struc", str(sid), )

        #########  OPT  #########
        if self.model.read.general.gfn2opt and len(ele_l) > 1:
            wd = create_folder(join(comdir, "xtb"))
            write_xyz(coordinate, ele_l, "coord")
            staus = CallxTB(wd, custom_settings={"charge": charge, "uhf": m - 1}).opt()
            if staus == 0:
                coordinate, ele_l = read_coord(join(wd, "xtbopt.xyz"))
            else:
                print(f'Warning: struc {sid} opt at xtb failed, continue with unopt struc')

        # do DFT opt if required
        forcedftopt = self.model.read.general.dftopt
        generalopt = self.model.orca_setting.dftopt
        if (forcedftopt or generalopt) and len(ele_l) > 1:
            wd = create_folder(join(comdir, "orca"))
            write_xyz(coordinate, ele_l, "coord")
            staus = CallORCA(wd, custom_settings={"charge": charge, "m": m, "caltyp": "SP OPT"}).call()
            if staus == 0:
                coordinate, ele_l = read_coord(join(wd, "orca.xyz"))
            else:
                print(f'Warning: struc {sid} opt at ORCA failed, continue with unopt struc')

        #########  SP  #########

        ## xtb sp
        if not exists(join(comdir, "xtb", "xtbopt.xyz")):
            wd = create_folder(join(comdir, "xtb"))
            write_xyz(coordinate, ele_l, "coord")
            staus = CallxTB(wd, custom_settings={"charge": charge, "uhf": m - 1}).sp()
            if staus != 0:
                print(f'Warning: struc {sid} sp at xtb failed, continue anyhow')

        ## DFT or gxtb sp
        # priority gxtb > orca
        if self.model.searching_setting.general.sp_method == "dft":
            if not (forcedftopt or generalopt) or len(ele_l) == 1:
                wd = create_folder(join(comdir, "orca"))
                write_xyz(coordinate, ele_l, "coord")
                staus = CallORCA(wd, custom_settings={"charge": charge, "m": m}).call()
                if staus != 0:
                    print(f'Warning: struc {sid} sp at ORCA failed, continue anyhow')
        else:
            # gxtb sp
            wd = create_folder(join(comdir, "gxtb"))
            write_xyz(coordinate, ele_l, "coord")
            par = CallgXTB(wd, custom_settings={"charge": charge, "uhf": m - 1})
            # write .charg and .uhf file
            with open(".CHRG", "w") as file:
                file.write(f"{charge}\n")
            with open(".UHF", "w") as file:
                file.write(f"{m - 1}\n")
            status = par.sp()

        # #########  TRV  #########
        # wd = create_folder(join(Setting.strucwd, strucid, "xtb", "trv"))
        # write_xyz(coord, "coord")
        # staus = CallxTB(wd, custom_settings={"charge": charge, "uhf": m - 1}).trv()
        # if staus != 0:
        #     print(f'Warning: struc {strucid} trv failed, continue anyhow')

        #########  SOLV  #########
        # wd = create_folder(join(Setting.strucwd, strucid, "xtb", "solv"))
        # write_xyz(coord, "coord")
        # staus = CallxTB(wd, custom_settings={"charge": charge, "uhf": m - 1}).solv()
        # if staus != 0:
        #     print(f'Warning: struc {strucid} solv failed, continue anyhow')

        # #########  get res  #########
        struc_info, complete = get_strucinfo_from_path(comdir)

        return struc_info

    def get(self, sid):
        return self.data[sid]


class Allrxn:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance.data = mydict({})
        return cls._instance

    def get_value(self, target_key, default_value=None):
        value = self.data.get_value(target_key, default_value)
        return value

    def set_value(self, target_key, new_value):
        self.data.set_value(target_key, new_value)

    def compare_value(self, key1, key2):
        return self.data.compare_value(key1, key2)

    def get_all_keys(self):
        return list(self.data.keys())

    def get_layer_keys(self, layer=0):
        # cycle0: 0_0, 0_1, 0_2....
        # cycle1: 0_0_0, 0_0_1,0_0_2......
        allkeys = self.get_all_keys()
        return [_ for _ in allkeys if len(_.split("_")) == layer + 1]


setting = None
allstruc = None
allrxn = Allrxn()


def initialize(setting_instance):
    """Initialize global setting and allstruc with the provided Setting instance."""
    global setting, allstruc
    setting = setting_instance.model  # Assign the passed Setting instance
    allstruc = Allstruc(setting)
