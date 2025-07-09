import sys
import toml
import threading
import os
from os.path import join, exists
from aRST.script.toolbox import create_folder, mydict, read_coord, get_strucinfo_from_path, write_xyz, get_all_keys
from aRST.script.callsys import CallxTB, CallORCA
from aRST.script.check_input import checking_parameter, checking_env


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
        # changing all str parameter to lower cases
        all_keys = get_all_keys(self.config)
        for key in all_keys:
            val = self.get_value(key)
            if isinstance(val, str) and not "geom" in key:
                self.set_value(key, val.lower())

        # checking up parameter format

        self.config = checking_parameter(self.config)
        checking_env(self.config)

        self.wd0 = self.get_value("searching_setting.general.wd", os.path.abspath(os.path.dirname(infile)))
        self.scanwd = join(self.wd0, "scan")
        self.strucwd = join(self.wd0, "struc")
        self.bufferwd = join(self.wd0, "buffer_struc")
        # clean up buffer struc
        if exists(self.bufferwd):
            import shutil
            shutil.rmtree(self.bufferwd)

    def get_value(self, target_key, default_value=None):
        value = mydict(self.config).get_value(target_key, default_value)
        return value

    def set_value(self, target_key, new_value):
        mydict(self.config).set_value(target_key, new_value)


class Allstruc:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, Setting=None):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance._initialize(Setting)
        return cls._instance

    def _initialize(self, Setting):
        self.number_of_struc = mydict(
            {})  # only used in nanoreactor for recording change of number of structures in sys
        self.number_of_struc["0"] = []
        # self.number_of_struc = {"layer":[{"strucid": strucid,
        #                                   "num": num}, ...]}
        self.data = mydict({})
        self.reactants_ini = []
        # read geom input
        ## reactants
        reactants = [_ for _ in Setting.get_value("geom").keys() if "reactant" in _]
        for id, struc in enumerate(reactants):
            coord = read_coord(join(Setting.wd0, Setting.get_value(f"geom.{struc}.name")))
            c = Setting.get_value(f"geom.{struc}.charge")
            m = Setting.get_value(f"geom.{struc}.multiplicity")
            if all(v is not None for v in [coord, c, m]):
                # set up new id anyhow
                self.data.set_value(f"{str(id)}.struc_info.general.coord", coord)
                self.data.set_value(f"{str(id)}.struc_info.general.charge", c)
                self.data.set_value(f"{str(id)}.struc_info.general.multiplicity", m)
            else:
                raise ValueError("Can't read reactants info from input, check if correct coord name, charge, m!")
            nstruc = Setting.get_value(f"geom.{struc}.number", 1)
            self.data.set_value(f"{str(id)}.struc_info.general.number", nstruc)
            self.number_of_struc["0"].append({"strucid": str(id), "num": int(nstruc)})
            for _ in range(nstruc):
                self.reactants_ini.append(str(id))

        ## ref products
        self.refproducts = []
        self.data_refp = mydict({})
        products = [_ for _ in Setting.get_value("geom").keys() if "product" in _]
        for id, struc in enumerate(products):
            coord = read_coord(join(Setting.wd0, Setting.get_value(f"geom.{struc}.name")))
            c = Setting.get_value(f"geom.{struc}.charge")
            m = Setting.get_value(f"geom.{struc}.multiplicity")
            if all(v is not None for v in [coord, c, m]):
                self.data_refp.set_value(f"{str(id)}.struc_info.general.coord", coord)
                self.data_refp.set_value(f"{str(id)}.struc_info.general.charge", c)
                self.data_refp.set_value(f"{str(id)}.struc_info.general.multiplicity", m)
            else:
                raise ValueError("Can't read refproducts info from input, check if correct coord name, charge, m!")
            for _ in range(Setting.get_value(f"geom.{struc}.number", 1)):
                self.refproducts.append(str(id))

        # collect struc refine info
        if Setting.get_value(f"read.general.jumprefine", False):
            # direct read struc refine information from foders 0, 1, etc
            for id, struc in enumerate(reactants):
                try:
                    struc_info = get_strucinfo_from_path(join(Setting.wd0, "struc", str(id)))
                    self.set_value_from_struc_info(str(id), struc_info)
                except:
                    print(f"{str(id)} strucinfo doesn't exist!, start to refine it")
                    coord = read_coord(join(Setting.wd0, Setting.get_value(f"geom.{struc}.name")))
                    c = Setting.get_value(f"geom.{struc}.charge")
                    m = Setting.get_value(f"geom.{struc}.multiplicity")
                    struc_info = self.refine_struc_info(coord, c, m, Setting, str(id))
                    self.set_value_from_struc_info(str(id), struc_info)


        else:
            for id, struc in enumerate(reactants):
                coord = read_coord(join(Setting.wd0, Setting.get_value(f"geom.{struc}.name")))
                c = Setting.get_value(f"geom.{struc}.charge")
                m = Setting.get_value(f"geom.{struc}.multiplicity")
                struc_info = self.refine_struc_info(coord, c, m, Setting, str(id))
                self.set_value_from_struc_info(str(id), struc_info)

    def set_value_from_struc_info(self, strucid, struc_info):
        # test if this reactant is already in data
        mb = struc_info["mb"]
        for searchid in list(self.data.keys()):
            if searchid != strucid:
                search_mb = self.data.get_value(f"{searchid}.struc_info.general.mb")
                if mb == search_mb:
                    del self.data[strucid]
                    self.reactants_ini.remove(strucid)
                    self.reactants_ini.append(searchid)
                    return
        # check if need to change current strucid
        for _ in range(int(strucid)):
            if str(_) not in list(self.data.keys()):
                del self.data[strucid]
                strucid = str(_)
                break

        self.data.set_value(f"{strucid}.struc_info.general.coord", struc_info["coord"])
        self.data.set_value(f"{strucid}.struc_info.general.nel", struc_info["nel"])
        self.data.set_value(f"{strucid}.struc_info.general.multiplicity", struc_info["m"])
        self.data.set_value(f"{strucid}.struc_info.general.charge", struc_info["c"])
        self.data.set_value(f"{strucid}.struc_info.general.atcharge", struc_info["atcharge"])
        self.data.set_value(f"{strucid}.struc_info.general.mb", struc_info["mb"])
        self.data.set_value(f"{strucid}.struc_info.general.atomic_connection",
                            (struc_info["constrain_l"], struc_info["ele_l"]))

        self.data.set_value(f"{strucid}.struc_info.energy.sp_gfn", struc_info["Esp_gfn"])
        self.data.set_value(f"{strucid}.struc_info.energy.sp_dft", struc_info["Esp_orca"])
        # self.data.set_value(f"{strucid}.struc_info.energy.trv", struc_info["Etrv_gfn"])
        # self.data.set_value(f"{strucid}.struc_info.energy.solv", struc_info["Esolv_gfn"])

    @staticmethod
    def refine_struc_info(coord, charge, m, Setting, strucid=0, wd0=None):

        #########  OPT  #########

        # do gfn opt if required
        if Setting.get_value("read.general.gfn2opt", False) and coord.shape[0] > 1:
            if not wd0:
                wd = create_folder(join(Setting.strucwd, strucid, "xtb"))
            else:
                wd = create_folder(join(wd0, "xtb"))
            write_xyz(coord, "coord")
            staus = CallxTB(wd, custom_settings={"charge": charge, "uhf": m - 1}).opt()
            if staus == 0:
                coord = read_coord(join(wd, "xtbopt.xyz"))
            else:
                print(f'Warning: struc {strucid} opt at xtb failed, continue with unopt struc')

        # do DFT opt if required
        forcedftopt = Setting.get_value("read.general.dftopt", False)
        if (forcedftopt or Setting.get_value("orca_setting.dftopt", False)) and coord.shape[0] > 1:
            if not wd0:
                wd = create_folder(join(Setting.strucwd, strucid, "orca"))
            else:
                wd = create_folder(join(wd0, "orca"))
            write_xyz(coord, "coord")
            staus = CallORCA(wd, custom_settings={"charge": charge, "m": m, "caltyp": "SP OPT"}).call()
            if staus == 0:
                coord = read_coord(join(wd, "orca.xyz"))
            else:
                print(f'Warning: struc {strucid} opt at ORCA failed, continue with unopt struc')

        #########  SP  #########

        ## xtb sp
        if not wd0:
            if not exists(join(Setting.strucwd, strucid, "xtb", "xtbopt.xyz")):
                wd = create_folder(join(Setting.strucwd, strucid, "xtb"))
                write_xyz(coord, "coord")
                staus = CallxTB(wd, custom_settings={"charge": charge, "uhf": m - 1}).sp()
                if staus != 0:
                    print(f'Warning: struc {strucid} sp at xtb failed, continue anyhow')
        else:
            wd = create_folder(join(wd0, "xtb"))
            write_xyz(coord, "coord")
            staus = CallxTB(wd, custom_settings={"charge": charge, "uhf": m - 1}).sp()
            if staus != 0:
                print(f'Warning: struc {strucid} sp at xtb failed, continue anyhow')

        ## orca sp
        sp_orca = None
        if not (forcedftopt or Setting.get_value("orca_setting.dftopt", False)) or coord.shape[0] == 1:
            if not wd0:
                wd = create_folder(join(Setting.strucwd, strucid, "orca"))
            else:
                wd = create_folder(join(wd0, "orca"))
            write_xyz(coord, "coord")
            staus = CallORCA(wd, custom_settings={"charge": charge, "m": m}).call()
            if staus != 0:
                print(f'Warning: struc {strucid} sp at ORCA failed, continue anyhow')

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

        #########  get res  #########
        if not wd0:
            struc_info = get_strucinfo_from_path(join(Setting.strucwd, strucid))

            return struc_info

    def get_value(self, target_key, default_value=None):
        value = self.data.get_value(target_key, default_value)
        return value

    def set_value(self, target_key, new_value):
        self.data.set_value(target_key, new_value)

    def compare_value(self, key1, key2):
        return self.data.compare_value(key1, key2)


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
    setting = setting_instance  # Assign the passed Setting instance
    allstruc = Allstruc(setting)
