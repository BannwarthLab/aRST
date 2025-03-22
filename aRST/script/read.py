import os
import shutil
import sys

from aRST.script.dictionaryInspector import DictionaryInspector
from aRST.script.toolbox import get_nelectron, create_folder, get_energy, write_xyz, read_topomol, get_molbarcode
from aRST.script.callsys import CallxTB, CallORCA
import pandas as pd

def read_geom_d(wd,setting_d):
    geom_d = setting_d["geom"]
    setting_inspertor = DictionaryInspector(setting_d)
    # read struc info and generate STRUC_rec
    STRUC_rec = {}
    REFSTRUC_rec = {}
    REACTION_rec = {}
    STRUC_inspector = DictionaryInspector(STRUC_rec)
    REFSTRUC_inspector = DictionaryInspector(REFSTRUC_rec)
    reactants_ini = []
    refproducts_ini = []
    rcount, pcount = 0, 0
    jumprefine = setting_inspertor.get_value(f"read.general.jumprefine", False)

    if jumprefine:
        nreactants = len([_ for _ in geom_d.keys() if "reactant" in _])
        for n in range(nreactants):
            struckn = list(geom_d.keys())[n]
            STRUC_rec = initial_struc_info(wd, geom_d, str(n), struckn, setting_d, STRUC_rec,
                                           jumprefine=jumprefine)
            reactants_ini.append(str(n))
        return reactants_ini,refproducts_ini,STRUC_rec,REFSTRUC_rec,REFSTRUC_rec,REACTION_rec
    for struckn in geom_d.keys():
        number = setting_inspertor.get_value(f"geom.{struckn}.number", 1)
        if "reactant" in struckn:
            if not reactants_ini:
                STRUC_rec = initial_struc_info(wd, geom_d, str(rcount), struckn, setting_d, STRUC_rec,
                                               jumprefine=jumprefine)
                reactants_ini = [str(rcount) for _ in range(setting_inspertor.get_value(f"geom.{struckn}.number", 1))]
                rcount += 1
            else:
                tmp_STRUC_rec = initial_struc_info(wd, geom_d, str(rcount), struckn, setting_d, {},
                                                   jumprefine=jumprefine)
                tmp_mb = DictionaryInspector(tmp_STRUC_rec).get_value(f"{rcount}.struc_info.general.mb")
                addtriger = True
                if tmp_mb != "unknow":
                    for struc in STRUC_rec.keys():
                        if STRUC_inspector.get_value(f"{struc}.struc_info.general.mb") == tmp_mb:
                            addtriger = False
                            for _ in range(setting_inspertor.get_value(f"geom.{struckn}.number", 1)):
                                reactants_ini.append(struc)
                if addtriger:
                    STRUC_rec[str(rcount)] = tmp_STRUC_rec[str(rcount)]
                    for _ in range(setting_inspertor.get_value(f"geom.{struckn}.number", 1)):
                        reactants_ini.append(str(rcount))
                    rcount += 1
                else:
                    # del the tmp calculation folder
                    tmp_struc_wd = os.path.join(wd, "struc", str(rcount), "xtb")
                    if os.path.exists(tmp_struc_wd):
                        shutil.rmtree(tmp_struc_wd, ignore_errors=True)
        elif "product" in struckn:
            if not refproducts_ini:
                REFSTRUC_rec = initial_struc_info(wd, geom_d, str(pcount), struckn, setting_d, REFSTRUC_rec,
                                                  jumprefine=jumprefine)
                refproducts_ini = [str(pcount) for _ in range(setting_inspertor.get_value(f"geom.{struckn}.number", 1))]
                pcount += 1
            else:
                tmp_REFSTRUC_rec = initial_struc_info(wd, geom_d, str(pcount), struckn, setting_d, {},
                                                      jumprefine=jumprefine)
                tmp_mb = DictionaryInspector(tmp_REFSTRUC_rec).get_value(f"{pcount}.struc_info.general.mb")
                addtriger = True
                if tmp_mb != "unknow":
                    for struc in REFSTRUC_rec.keys():
                        if REFSTRUC_inspector.get_value(f"{struc}.struc_info.general.mb") == tmp_mb:
                            addtriger = False
                            for _ in range(setting_inspertor.get_value(f"geom.{struckn}.number", 1)):
                                refproducts_ini.append(struc)
                if addtriger:
                    REFSTRUC_rec[str(pcount)] = tmp_REFSTRUC_rec[str(rcount)]
                    for _ in range(setting_inspertor.get_value(f"geom.{struckn}.number", 1)):
                        refproducts_ini.append(str(pcount))
                    pcount += 1
    return reactants_ini,refproducts_ini,STRUC_rec,REFSTRUC_rec,REFSTRUC_rec,REACTION_rec


def initial_struc_info(wd0,geom_d, strucid, keyname,setting_d,struc_rec = {},jumprefine=False):
    os.chdir(wd0)
    geom_inspector = DictionaryInspector(geom_d)
    STRUC_inspector = DictionaryInspector(struc_rec)
    coord_fn = geom_inspector.get_value(f"{keyname}.name")
    c, m = geom_inspector.get_value(f"{keyname}.charge"), geom_inspector.get_value(f"{keyname}.multiplicity")
    if all(v is not None for v in [coord_fn, c , m]):
        STRUC_inspector.set_value(f"{strucid}.struc_info.general.coord",
                                  pd.read_table(coord_fn, header=None, skiprows=2, sep='\s+'))
        STRUC_inspector.set_value(f"{strucid}.struc_info.general.charge", c)
        STRUC_inspector.set_value(f"{strucid}.struc_info.general.multiplicity", m)

        if geom_inspector.get_value(f"{keyname}.at") is not None:
            STRUC_inspector.set_value(f"{strucid}.struc_info.general.at", geom_inspector.get_value(f"{keyname}.at"))
        #
        if not jumprefine:
            struc_rec = refine_struc_info(wd0,strucid,struc_rec,setting_d)
        return struc_rec
    else:
        raise Exception("Can't read struc info from input file.")


def refine_struc_info(wd0,strucid, struc_rec, setting_d):
    struc_inspector = DictionaryInspector(struc_rec)
    setting_inspector = DictionaryInspector(setting_d)
    # get bisical properties
    coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
    charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
    m = struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity")

    # do QC calculations to get more struc properties
    nel, sp_xtb = None,None
    constrain_l, ele = None,None

    # do gfn opt if required
    if setting_inspector.get_value("read.general.gfn2opt", False) and coord.shape[0] > 1:
        wd = create_folder(os.path.join(wd0, "struc", strucid, "xtb"))
        write_xyz(coord, "coord")

        staus = CallxTB(wd, {"method": setting_inspector.get_value("xtb_setting.method", "gfn2"),
                             "coordfn": "coord",
                             "optlevel": setting_inspector.get_value("xtb_setting.optlevel", ""),
                             "charge": charge,
                             "uhf": m - 1,
                             "alpb": setting_inspector.get_value("xtb_setting.alpb", "h2o"),
                             "etemp": setting_inspector.get_value("xtb_setting.struc_etemp", "4000")
                             }).opt()
        if staus == 0:
            coord = pd.read_table(os.path.join(wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
            nel = get_nelectron(os.path.join(wd, 'xtb.out'))
            sp_xtb = get_energy(os.path.join(wd, 'ENERGY'))
            if os.path.exists(os.path.join(wd, 'xtbtopo.mol')):
                constrain_l, ele = read_topomol(os.path.join(wd, 'xtbtopo.mol'))
            else:
                constrain_l, ele = [], []
        else:
            nel, sp_xtb = None, None
            constrain_l, ele = [], [str(coord.iloc[0, 0])]
            print(f'struc {strucid} opt at {setting_inspector.get_value("xtb_setting.method", "gfn2")} failed')

    # do DFT opt if required
    forcedftopt = setting_inspector.get_value("read.general.dftopt", False)
    if (forcedftopt or setting_inspector.get_value("orca_setting.dftopt") )and coord.shape[0] > 1:
        dftopt = True
        # do dft opt
        wd = create_folder(os.path.join(wd0,"struc", strucid, "orca"))
        write_xyz(coord, "coord")
        staus = CallORCA(wd, {"functional": setting_inspector.get_value("orca_setting.functional", "r2scan-3c"),
                              "basis": setting_inspector.get_value("orca_setting.basis", " "),
                              "charge": charge,
                              "m": m,
                              "fn": "coord",
                              "caltyp": "SP OPT",
                              "caltyp2": "RKS" if m == 1 else "UKS",
                              "solmode": setting_inspector.get_value("orca_setting.solmode", ""),
                              "ncpu": setting_inspector.get_value("orca_setting.ncpu", 1)}).call()
        if staus == 0:
            coord = pd.read_table(os.path.join(wd, 'orca.xyz'), header=None, skiprows=2, sep='\s+')

        else:
            print(f'struc {strucid} opt at '
                            f'{setting_inspector.get_value("orca_setting.functional", "r2scan-3c")} failed')

    ### do sp
    # do xtb sp
    if not all([_ is not None for _ in [nel,sp_xtb,ele,constrain_l]]):
        wd = create_folder(os.path.join(wd0, "struc", strucid, "xtb"))
        write_xyz(coord, "coord")
        staus = CallxTB(wd, {"method": setting_inspector.get_value("xtb_setting.method", "gfn2"),
                             "coordfn": "coord",
                             "optlevel": setting_inspector.get_value("xtb_setting.optlevel", ""),
                             "charge": charge,
                             "uhf": m - 1,
                             "alpb": setting_inspector.get_value("xtb_setting.alpb", "h2o"),
                             "etemp": setting_inspector.get_value("xtb_setting.struc_etemp", "4000")
                             }).sp()
        if staus == 0:
            nel = get_nelectron(os.path.join(wd, 'xtb.out'))
            sp_xtb = get_energy(os.path.join(wd, 'ENERGY'))
            if os.path.exists(os.path.join(wd, 'xtbtopo.mol')):
                constrain_l, ele = read_topomol(os.path.join(wd, 'xtbtopo.mol'))
            elif coord.shape[0]==1:
                constrain_l=[]
                ele = [coord.iloc[0,0]]
            else:
                constrain_l, ele = [], []
        else:
            nel, sp_xtb = None, None
            raise Exception(
                f"struc {strucid} sp at {setting_inspector.get_value('xtb_setting.method', 'gfn2')} failed")


    # # do orca sp
    sp_orca=0
    if not setting_inspector.get_value("orca_setting.jumpsp", False):
        wd = create_folder(os.path.join(wd0, "struc", strucid, "orca"))
        write_xyz(coord, "coord")
        staus = CallORCA(wd, {"functional": setting_inspector.get_value("orca_setting.functional", "r2scan-3c"),
                              "basis": setting_inspector.get_value("orca_setting.basis", " "),
                              "charge": charge,
                              "m": m,
                              "fn": "coord",
                              "caltyp": "SP",
                              "caltyp2": "RKS" if m == 1 else "UKS",
                              "solmode": setting_inspector.get_value("orca_setting.solmode", ""),
                              "ncpu": setting_inspector.get_value("orca_setting.ncpu", 1)}).call()
        if staus == 0:
            sp_orca = get_energy(os.path.join(wd, 'ENERGY'))
        else:
            sp_orca = None
            raise Exception(f'struc {strucid} sp at '
                            f'{setting_inspector.get_value("orca_setting.functional", "r2scan-3c")} failed')
    endxyz = coord

    # getting molbar code
    mb = get_molbarcode(endxyz, setting_inspector.get_value("searching_setting.general.mb_levlel"))

    # getting trv and solv
    if coord.shape[0]>1:
        wd = create_folder(os.path.join(wd0,"struc", strucid, "xtb", "trv"))
        write_xyz(coord, "coord")
        staus = CallxTB(wd, {"method": setting_inspector.get_value("xtb_setting.method", "gfn2"),
                             "coordfn": "coord",
                             "optlevel": setting_inspector.get_value("xtb_setting.optlevel", ""),
                             "charge": charge,
                             "uhf": m - 1,
                             "alpb": setting_inspector.get_value("xtb_setting.alpb", "h2o"),
                             "etemp": setting_inspector.get_value("xtb_setting.struc_etemp", "4000")
                             }).trv()
        etrv = get_energy(os.path.join(wd, 'GTRV')) if staus == 0 else 0
    else:etrv=0


    wd = create_folder(os.path.join(wd0,"struc", strucid, "xtb", "solv"))
    write_xyz(coord, "coord")
    staus = CallxTB(wd, {"method": setting_inspector.get_value("xtb_setting.method", "gfn2"),
                         "coordfn": "coord",
                         "optlevel": setting_inspector.get_value("xtb_setting.optlevel", ""),
                         "charge": charge,
                         "uhf": m - 1,
                         "alpb": setting_inspector.get_value("xtb_setting.alpb", "h2o"),
                         "etemp": setting_inspector.get_value("xtb_setting.struc_etemp", "4000")
                         }).solv()
    esolv = get_energy(os.path.join(wd, 'GSOLV')) if staus == 0 else 0

    # collect data
    struc_inspector.set_value(f"{strucid}.struc_info.general.coord", coord)
    struc_inspector.set_value(f"{strucid}.struc_info.general.nel", nel)
    struc_inspector.set_value(f"{strucid}.struc_info.general.mb", mb)
    struc_inspector.set_value(f"{strucid}.struc_info.general.atomic_connection", (constrain_l, ele))
    struc_inspector.set_value(f"{strucid}.struc_info.energy.sp_gfn", sp_xtb)
    struc_inspector.set_value(f"{strucid}.struc_info.energy.sp_dft", sp_orca)
    struc_inspector.set_value(f"{strucid}.struc_info.energy.trv", etrv)
    struc_inspector.set_value(f"{strucid}.struc_info.energy.solv", esolv)
    # struc_inspector.set_value(f"{strucid}.struc_info.energy.G_gfn", sp_xtb + etrv + esolv)
    struc_inspector.set_value(f"{strucid}.struc_info.energy.G_dft", sp_orca + etrv + esolv)
    return struc_rec
