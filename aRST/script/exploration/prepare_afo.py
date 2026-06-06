import os
import pandas as pd
import numpy as np
from aRST.script.toolbox import create_folder, write_xyz, read_bo, read_afo_csv, timeit
from os.path import join, exists, dirname
from aRST.script.callsys import CallxTB, CallORCA, CallTerachem, CallsTDA, CallMultiwfn, CallgXTB
import aRST.script.initialization as storage
import shutil


class ReactivityPreparer:
    def __init__(self, sid_list):
        self.sid_list = sid_list

    @staticmethod
    def check_featrues(sid):
        # check if coord, nel, afo, mbo exist
        s = storage.allstruc.get(sid)
        coordinate = s.general.coordinate
        nel = s.general.nel
        c = s.general.charge
        if not all(v is not None for v in [coordinate, nel, c]):
            print(
                f"Warning: {sid} doesn't have enough coord, nel, c, this struc will be neglet in the reaction searching")

        afo = s.reactivity.afo
        if afo is None:
            return False
        else:
            return True

    @timeit("prepare AFO features")
    def ensure_AFO_ready(self, afosubwd=None):
        self.prepare_AFO(afosubwd=afosubwd)

    def ensure_weakH_ready(self):
        self.prepare_weakH_list()

    def prepare_AFO(self, afosubwd=None, force=False):
        for sid in self.sid_list:
            if not force and self.check_featrues(sid):
                continue

            s_generalentry = storage.allstruc.get(sid).general
            s_reactivityentry = storage.allstruc.get(sid).reactivity
            c = s_generalentry.charge
            ele_l = s_generalentry.ele_l

            if len(ele_l) == 1 and ele_l[0].strip().upper() == "H":
                mbo = None
                afo = pd.DataFrame(
                    columns=['index', 'HOAO_a (eV)', 'HOAO_b (eV)', 'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom'])
                if c == 1:
                    # H+
                    afo.loc[0, :] = [0, -1e5, -1e5, -1e5, -1e5, 1]

                elif c == -1:
                    # H-
                    afo.loc[0, :] = [0, 1e5, 1e5, 1e5, 1e5, 1]
                else:
                    # H radical
                    afo.loc[0, :] = [0, 1e5, -1e5, 1e5, -1e5, 1]
                s_reactivityentry.afo = afo
                s_reactivityentry.bo = mbo
                continue

            if str.lower(storage.setting.searching_setting.reactivity.afo_method) == "gfn2":
                afowd, status = getAFO(strucid=sid).get_from_xtb(afosubwd)
            elif str.lower(storage.setting.searching_setting.reactivity.afo_method) == "terachem":
                afowd, status = getAFO(strucid=sid).get_from_tc(afosubwd)
            else:
                # defalut as orca sp
                afowd, status = getAFO(strucid=sid).get_from_orca(afosubwd)

            # read afo
            afo = read_afo_csv(join(afowd, "ml_feature.csv"))

            # ###### del H radical beta afo features (no beta electron at all) ######
            # if nel == 1:
            #     afo.loc[:, "HOAO_b (eV)"] = -np.inf
            #     afo.loc[:, "LUAO_b (eV)"] = np.inf

            s_reactivityentry.afo = afo

            # set bo
            for file in [join(afowd, "wboe"), join(afowd, "bde"), join(afowd, "embo")]:
                mbo = read_bo(file)
                if mbo is not None:
                    break
            s_reactivityentry.bo = mbo

    def prepare_weakH_list(self):
        for sid in self.sid_list:
            s = storage.allstruc.get(sid)
            weakH = s.reactivity.weakH
            if weakH is not None:
                continue

            ele_l = s.general.ele_l
            allH_idx = [i for i, ele in enumerate(ele_l) if ele.lower() == "h"]

            if len(allH_idx)==len(ele_l): # only 1 H atom
                s.reactivity.weakH = allH_idx
                continue

            topo_l = s.general.atomic_connection
            weakHl = []

            for Hidx in allH_idx:
                Hidx_1 = Hidx + 1
                mbo = next(
                    (
                        item["mbo"]
                        for item in topo_l
                        if Hidx_1 in item["atoms"]
                    ),
                    None
                )

                if mbo is not None and mbo < 0.9:
                    weakHl.append(Hidx)


            s.reactivity.weakH = weakHl
        return


class getAFO:
    def __init__(self, strucid=0, coordinate=None, ele_l=None
                 , c=None, m=None):
        if coordinate is None:
            s = storage.allstruc.get(strucid)
            self.coordinate = s.general.coordinate
            self.ele_l = s.general.ele_l
            self.c = s.general.charge
            self.m = s.general.multiplicity
        else:
            self.coordinate = coordinate
            self.ele_l = ele_l
            self.c = c
            self.m = m
        self.struc = str(strucid)
        self.strucwd = storage.setting.searching_setting.general.strucwd

    def get_from_xtb(self, subfolder=None):
        afowd = create_folder(join(self.strucwd, self.struc, "xtb", "afo")) \
            if not subfolder else create_folder(join(self.strucwd, self.struc, subfolder, "xtb", "afo"))
        rc = 0
        if not exists(join(afowd, "ml_feature.csv")):
            if 'XTBML_BIN' in os.environ and os.environ['XTBML_BIN']:
                # modified xtb exiest
                write_xyz(self.coordinate, self.ele_l, "coord", c=self.c, m=self.m)
                rc = CallxTB(afowd, custom_settings={"charge": self.c, "uhf": self.m - 1},mlversion=True).afo()
            else:
                # use normal xtb version to generate molden file
                write_xyz(self.coordinate, self.ele_l, "coord", c=self.c, m=self.m)
                rc = CallxTB(afowd, custom_settings={"charge": self.c, "uhf": self.m - 1}).molden()
                if rc == 0:
                    rc = CallsTDA(afowd, "molden.input").callml(0)
                else:
                    rc = 1
        return afowd, rc

    def get_from_gxtb_molden(self, subfolder=None):
        afowd = create_folder(join(self.strucwd, self.struc, "xtb", "afo")) \
            if not subfolder else create_folder(join(self.strucwd, self.struc, subfolder, "xtb", "afo"))
        rc = 0
        if not exists(join(afowd, "ml_feature.csv")):
            write_xyz(self.coordinate, self.ele_l, "coord", c=self.c, m=self.m)
            rc = CallgXTB(afowd, custom_settings={"charge": self.c, "uhf": self.m - 1}).afo_molden()
            if rc == 0:
                rc = CallsTDA(afowd, "molden.input").callml(0)
            else:
                rc = 1

        return afowd, rc

    def get_from_tc(self, subfolder=None):
        afowd = join(self.strucwd, self.struc, "stda") \
            if not subfolder else create_folder(join(self.strucwd, self.struc, subfolder, "stda"))

        rc = 0
        if not exists(join(afowd, "ml_feature.csv")):
            tcwd = create_folder(join(self.strucwd, self.struc, "tc"))
            write_xyz(self.coordinate, self.ele_l, "coord", c=self.c, m=self.m)
            inputcmd_d = {"coordfn": "coord",
                          "c": self.c,
                          "m": self.m,
                          "functional": storage.setting.searching_setting.reactivity.afo_tc_functional,
                          "basis": storage.setting.searching_setting.reactivity.afo_tc_basis,
                          "epsilon": storage.setting.searching_setting.reactivity.afo_tc_epsilon,
                          }

            rc = CallTerachem(wd=tcwd, command_d=inputcmd_d).run()

            if rc == 0:
                afowd = create_folder(join(self.strucwd, self.struc, "stda")) \
                    if not subfolder else create_folder(join(self.strucwd, self.struc, subfolder, "stda"))
                shutil.copyfile(join(tcwd, "scr.coord", "coord.molden"),
                                join(afowd, "coord.molden"))
                rc = CallsTDA(afowd, "coord.molden").callml(3)
        return afowd, rc

    def get_from_orca(self, subfolder=None):
        afowd = create_folder(join(self.strucwd, self.struc, "stda")) \
            if not subfolder else create_folder(join(self.strucwd, self.struc, subfolder, "stda"))
        rc = 0
        if not exists(join(afowd, "ml_feature.csv")):
            orcawd = create_folder(
                join(self.strucwd, self.struc, "orcaml")) if not subfolder else create_folder(
                join(self.strucwd, self.struc, subfolder, "orcaml"))
            write_xyz(self.coordinate, self.ele_l, "coord", c=self.c, m=self.m)
            staus = CallORCA(orcawd,
                             custom_settings={"charge": self.c, "m": self.m,
                                              "functional": storage.setting.searching_setting.reactivity.afo_orca_functional,
                                              "basis": storage.setting.searching_setting.reactivity.afo_orca_basis,
                                              "cpcm": storage.setting.searching_setting.reactivity.afo_orca_cpcm, }).call(
                type="sp",
                removegbw=False)

            import subprocess
            orca_2mkl_path = join(dirname(os.environ['ORCA_BIN']), "orca_2mkl")
            p = subprocess.Popen(f"{orca_2mkl_path} orca -molden", shell=True,
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL)
            p.wait()

            # transfom MOs from spherical coordinates to cartesian
            rc = CallMultiwfn(orcawd, "orca.molden.input").write_cartesian_molden()
            if rc == 0:
                self._safe_remove(join(orcawd, "orca.molden.input"))
                self._safe_remove(join(orcawd, "orca.gbw"))

                shutil.copyfile(join(orcawd, "orca.molden"),
                                join(afowd, "coord.molden"))
                self._safe_remove(join(orcawd, "orca.molden"))
                os.chdir(afowd)
                rc = CallsTDA(afowd, "coord.molden").callml(1)

                if rc == 0:
                    os.remove(join(afowd, "coord.molden"))
        return afowd, rc

    @staticmethod
    def _safe_remove(path):
        if exists(path):
            os.remove(path)
