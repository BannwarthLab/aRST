import os
from aRST.script.toolbox import create_folder, write_xyz
from os.path import join, exists, dirname
from aRST.script.callsys import CallxTB, CallORCA, CallTerachem, CallsTDA, CallMultiwfn,CallgXTB
import aRST.script.initialization as storage
import shutil


class getAFO:
    def __init__(self, strucid=0, coordinate=None, ele_l=None
                 , c=None, m=None):
        if coordinate is None:
            s = storage.allstruc.get(strucid)
            self.coordinate = s.struc_info.general.coordinate
            self.ele_l = s.struc_info.general.ele_l
            self.c = s.struc_info.general.charge
            self.m = s.struc_info.general.multiplicity
        else:
            self.coordinate = coordinate
            self.ele_l = ele_l
            self.c = c
            self.m = m
        self.struc = str(strucid)
        self.strucwd=storage.setting.searching_setting.general.strucwd

    def get_from_xtb(self, subfolder=None):
        afowd = create_folder(join(self.strucwd, self.struc, "xtb", "afo")) \
            if not subfolder else create_folder(join(self.strucwd, self.struc, subfolder, "xtb", "afo"))
        rc = 0
        if not exists(join(afowd, "ml_feature.csv")):
            if 'XTBML_BIN' in os.environ and os.environ['XTBML_BIN']:
                # modified xtb exiest
                write_xyz(self.coordinate, self.ele_l, "coord")
                rc = CallxTB(afowd, custom_settings={"charge": self.c, "uhf": self.m - 1}).afo()
            else:
                # use normal xtb version to generate molden file
                write_xyz(self.coordinate, self.ele_l, "coord")
                rc = CallxTB(afowd, custom_settings={"charge": self.c, "uhf": self.m - 1}).afo_molden()
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
            write_xyz(self.coordinate, self.ele_l, "coord")
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
            write_xyz(self.coordinate, self.ele_l, "coord")
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
            write_xyz(self.coordinate, self.ele_l, "coord")
            staus = CallORCA(orcawd,
                             custom_settings={"charge": self.c, "m": self.m,
                                              "functional": storage.setting.searching_setting.reactivity.afo_orca_functional,
                                              "basis": storage.setting.searching_setting.reactivity.afo_orca_basis,
                                              "solmode": storage.setting.searching_setting.reactivity.afo_orca_cpcm, }).call(
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

                if rc == 0:
                    os.remove(join(afowd, "coord.molden"))
        return afowd, rc
