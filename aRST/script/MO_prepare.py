import os

from aRST.script.dictionaryInspector import DictionaryInspector
from aRST.script.toolbox import create_folder, write_xyz
from os.path import join,exists
from aRST.script.callsys import CallxTB, CallORCA,CallTerachem,CallsTDA,CallMultiwfn
import shutil


"""
run sp calculations to get MOs and then get reactive features
"""
def _get_reactive_features(wd,struc,STRUC_rec,setting_d):
    STRUC_inspector = DictionaryInspector(STRUC_rec)
    setting_inspector = DictionaryInspector(setting_d)

    keyword1 = str.upper(setting_inspector.get_value("searching_setting.reactivity.ml_method"))

    if keyword1 == "GFN2": #
        mlwd = create_folder(join(wd, "struc", struc, "xtb", "ml"))
        if not exists(join(mlwd, "ml_feature.csv")):
            # call ml
            coord = STRUC_inspector.get_value(f"{struc}.struc_info.general.coord")
            charge = STRUC_inspector.get_value(f"{struc}.struc_info.general.charge")
            m = STRUC_inspector.get_value(f"{struc}.struc_info.general.multiplicity")
            write_xyz(coord, "coord")
            rc = CallxTB(mlwd, {"method": setting_d["xtb_setting"].get("method", "gfn2"),
                                "coordfn": "coord",
                                "optlevel": setting_d["xtb_setting"].get("optlevel", None),
                                "charge": charge,
                                "uhf": m - 1,
                                "alpb": setting_d["xtb_setting"].get("alpb", "h2o"),
                                "etemp": setting_d["xtb_setting"].get("struc_etemp", "4000"),
                                "xtbwd": setting_d["xtb_setting"].get("xtb_wd")
                                }).ml()
            if rc == 0:
                return mlwd
            else:
                raise Exception(f"struc {struc} xtb ml failed")
    else:
        # DFT based
        keyword2 = str.upper(setting_inspector.get_value("searching_setting.reactivity.MO_method"))
        if keyword2=="TERACHEM":
            # run sp in terachem
            tcwd = create_folder(join(wd, "struc", struc, "tc"))
            coord = STRUC_inspector.get_value(f"{struc}.struc_info.general.coord")
            write_xyz(coord, "coord")
            inputcmd_d = {"coordfn": "coord",
                          "c": STRUC_inspector.get_value(f"{struc}.struc_info.general.charge"),
                          "m": STRUC_inspector.get_value(f"{struc}.struc_info.general.multiplicity"),
                          "functional": setting_inspector.get_value("searching_setting.reactivity.ml_tc_functional"),
                          "basis": setting_inspector.get_value("searching_setting.reactivity.ml_tc_basis"),
                          "epsilon": setting_inspector.get_value("searching_setting.reactivity.ml_tc_epsilon"),
                          }

            rc = CallTerachem(wd=tcwd, command_d=inputcmd_d).run()

            if rc == 0:
                mlwd = create_folder(join(wd, "struc", struc, "stda"))
                shutil.copyfile(join(tcwd, "scr.coord", "coord.molden"),
                                join(mlwd, "coord.molden"))
                rc2 = CallsTDA(mlwd,"coord.molden").callml(setting_inspector.get_value("stda_setting.stda_wd"), 3)
                if rc2 == 0:
                    return mlwd
                else:
                    raise Exception(f"struc {struc} stda ml failed")

            else:
                raise Exception(f"struc {struc} tc molden failed")
        else:
            # run sp in ORCA
            orcawd = create_folder(join(wd, "struc", struc, "orcaml"))
            if not exists(join(orcawd,"orca.gbw")) and not exists(join(orcawd,"orca.molden")):
                coord = STRUC_inspector.get_value(f"{struc}.struc_info.general.coord")
                write_xyz(coord, "coord")
                m=STRUC_inspector.get_value(f"{struc}.struc_info.general.multiplicity")
                staus = CallORCA(orcawd, {"functional": setting_inspector.get_value("searching_setting.reactivity.ml_orca_functional"),
                                      "basis": setting_inspector.get_value("searching_setting.reactivity.ml_orca_basis"),
                                      "charge": STRUC_inspector.get_value(f"{struc}.struc_info.general.charge"),
                                      "m": m,
                                      "fn": "coord",
                                      "caltyp": "SP",
                                      "caltyp2": "RKS" if m == 1 else "UKS",
                                      "solmode": setting_inspector.get_value("searching_setting.reactivity.ml_orca_cpcm", ""),
                                      "ncpu": setting_inspector.get_value("orca_setting.ncpu", 1)}).call(removegbw=False)
                if staus==0:
                    pass
                else:raise Exception(f"struc {struc} orca ml failed")

                # get spherical coordinates type MOs
                import subprocess
                p = subprocess.Popen("orca_2mkl orca -molden", shell=True)
                p.wait()
                ### del .gbw file to keep i/o smaller
                os.remove(join(orcawd,"orca.gbw"))
                # transfom MOs from spherical coordinates to cartesian
                rc = CallMultiwfn(orcawd,"orca.molden.input").write_cartesian_molden(
                    setting_inspector.get_value("mwfn_setting.mwfn_wd"))
                if rc == 0:
                    os.remove(join(orcawd, "orca.molden.input"))
                    pass
                else:
                    raise Exception(f"struc {struc} Multiwfn transform failed")

            mlwd = create_folder(join(wd, "struc", struc, "stda"))
            shutil.copyfile(join(orcawd, "orca.molden"),
                            join(mlwd, "coord.molden"))
            rc2 = CallsTDA(mlwd, "coord.molden").callml(setting_inspector.get_value("stda_setting.stda_wd"), 1)
            if rc2 == 0:
                return mlwd
            else:
                raise Exception(f"struc {struc} stda ml failed")


    # elif keyword in ["gfn1,gfnff"]:  # need to call molden file first
    #     mlwd = create_folder(os.path.join(self.wd, "struc", struc, "stda"))
    #     if not os.path.exists(os.path.join(mlwd, "ml_feature.csv")):
    #         coord = STRUC_inspector.get_value(f"{struc}.struc_info.general.coord")
    #         charge = STRUC_inspector.get_value(f"{struc}.struc_info.general.charge")
    #         m = STRUC_inspector.get_value(f"{struc}.struc_info.general.multiplicity")
    #         write_xyz(coord, "coord")
    #         rc = CallxTB(mlwd, {"method": self.setting_d["xtb_setting"].get("method", "gfn2"),
    #                             "coordfn": "coord",
    #                             "optlevel": self.setting_d["xtb_setting"].get("optlevel", None),
    #                             "charge": charge,
    #                             "uhf": m - 1,
    #                             "alpb": self.setting_d["xtb_setting"].get("alpb", "h2o"),
    #                             "etemp": self.setting_d["xtb_setting"].get("struc_etemp", "4000")
    #                             }).get_moldeninput()
    #         if rc == 0:
    #             rc2 = CallsTDA("molden.input").callml(setting_inspector.get_value("stda_setting.stda_wd"), 0)
    #             if rc2 == 0:
    #                 pass
    #             else:
    #                 raise Exception(f"struc {struc} stda ml failed")
    #         else:
    #             raise Exception(f"struc {struc} xtb get molden file failed")
