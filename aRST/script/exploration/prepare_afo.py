import os
from aRST.script.toolbox import create_folder, write_xyz
from os.path import join, exists,dirname
from aRST.script.callsys import CallxTB, CallORCA, CallTerachem, CallsTDA, CallMultiwfn
import aRST.script.initialization as storage
import shutil


class getAFO:
    def __init__(self, strucid=0, setting=None, coord=None, c=None, m=None):
        if coord is None:
            self.coord = storage.allstruc.get_value(f"{strucid}.struc_info.general.coord")
            self.c = storage.allstruc.get_value(f"{strucid}.struc_info.general.charge")
            self.m = storage.allstruc.get_value(f"{strucid}.struc_info.general.multiplicity")
            self.setting = storage.setting

        else:
            self.coord = coord
            self.c = c
            self.m = m
            self.setting = setting
        self.struc = strucid

    def get_from_xtb(self,subfolder=None):
        afowd = create_folder(join(self.setting.strucwd, self.struc, "xtb", "afo")) \
            if not subfolder else create_folder(join(self.setting.strucwd, self.struc,subfolder, "xtb", "afo"))
        rc = 0
        if not exists(join(afowd, "ml_feature.csv")):
            write_xyz(self.coord, "coord")
            rc = CallxTB(afowd, custom_settings={"charge": self.c, "uhf": self.m - 1}).afo()
        return afowd, rc

    def get_from_tc(self,subfolder=None):
        afowd = join(self.setting.strucwd, self.struc, "stda") \
            if not subfolder else create_folder(join(self.setting.strucwd, self.struc,subfolder, "stda"))

        rc = 0
        if not exists(join(afowd, "ml_feature.csv")):
            tcwd = create_folder(join(self.setting.strucwd, self.struc, "tc"))
            write_xyz(self.coord, "coord")
            inputcmd_d = {"coordfn": "coord",
                          "c": self.c,
                          "m": self.m,
                          "functional": self.setting.get_value("searching_setting.reactivity.afo_tc_functional"),
                          "basis": self.setting.get_value("searching_setting.reactivity.afo_tc_basis"),
                          "epsilon": self.setting.get_value("searching_setting.reactivity.afo_tc_epsilon"),
                          }

            rc = CallTerachem(wd=tcwd, command_d=inputcmd_d).run()

            if rc == 0:
                afowd = create_folder(join(self.setting.strucwd,self.struc, "stda")) \
                    if not subfolder else create_folder(join(self.setting.strucwd, self.struc,subfolder, "stda"))
                shutil.copyfile(join(tcwd, "scr.coord", "coord.molden"),
                                join(afowd, "coord.molden"))
                rc = CallsTDA(afowd, "coord.molden").callml(3)
        return afowd, rc

    def get_from_orca(self,subfolder=None):
        afowd = create_folder(join(self.setting.strucwd, self.struc, "stda")) \
            if not subfolder else create_folder(join(self.setting.strucwd, self.struc,subfolder, "stda"))
        rc = 0
        if not exists(join(afowd, "ml_feature.csv")):
            orcawd = create_folder(join(self.setting.strucwd, self.struc, "orcaml")) if not subfolder else create_folder(join(self.setting.strucwd, self.struc,subfolder, "orcaml"))
            write_xyz(self.coord, "coord")
            staus = CallORCA(orcawd,
                             custom_settings={"charge": self.c, "m": self.m,
                                              "functional": self.setting.get_value(
                                                  "searching_setting.reactivity.afo_orca_functional"),
                                              "basis": self.setting.get_value(
                                                  "searching_setting.reactivity.afo_orca_basis"),
                                              "solmode": self.setting.get_value(
                                                  "searching_setting.reactivity.afo_orca_cpcm", ""), }).call(
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
                os.remove(join(orcawd, "orca.molden.input"))
                os.remove(join(orcawd, "orca.gbw"))

                shutil.copyfile(join(orcawd, "orca.molden"),
                                join(afowd, "coord.molden"))
                os.remove(join(orcawd, "orca.molden"))
                os.chdir(afowd)
                rc = CallsTDA(afowd, "coord.molden").callml(1)

                if rc ==0 and not storage.setting.get_value("searching_setting.general.keepmolden", False):
                    os.remove(join(afowd, "coord.molden"))
        return afowd, rc