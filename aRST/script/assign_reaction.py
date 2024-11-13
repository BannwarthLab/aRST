import os
import shutil

import numpy as np
import pandas as pd
from aRST.script.toolbox import create_folder, write_xyz, get_energy, seperate_fragment, read_topomol, get_molbarcode
from aRST.script.MO_prepare import _get_reactive_features
from aRST.script.dictionaryInspector import DictionaryInspector
from aRST.script.callsys import CallxTB, CallORCA


def ReactionSearch(wd, STRUC_rec, setting_d, REACTION_rec, step=0):
    setting_inspector = DictionaryInspector(setting_d)
    ncycle = setting_inspector.get_value("searching_setting.general.cycle", 1)
    print("ncycle=", ncycle)

    while step < ncycle:
        # find reaction pathes to be explored in this cycle
        if not REACTION_rec:
            print("begin cycle 0")
            # cycle 0
            reactants_l = [i for i in STRUC_rec.keys()]
            singlecycle = SingleCycleSearch(STRUC_rec, REACTION_rec, setting_d, wd)
            explore_dl = singlecycle.explore(reactants_l, "0")
            for reactiontype_d in explore_dl:
                if reactiontype_d:
                    for count, info_d in reactiontype_d.items():
                        # get products geom
                        info_d = singlecycle.follow_up_treatment_gfn(info_d)
                        # refine products info
                        info_d = singlecycle.follow_up_treatment_energy(info_d)
                        # now evaluate this reaction path
                        reaction_id = "0_" + str(count)
                        # update rec
                        STRUC_rec, REACTION_rec = singlecycle.evalueate_path(info_d, reaction_id)
                        singlecycle = SingleCycleSearch(STRUC_rec, REACTION_rec, setting_d, wd)

            step += 1
        else:
            print("begin cycle ", step)
            # cycle > 0
            node = [reaction_id for reaction_id in REACTION_rec.keys() if len(reaction_id.split("_")) == step + 1]
            print("node = ", node)
            for prefix_rpid in node:
                status = REACTION_rec[prefix_rpid][-1]
                if status:  # this reaction path is allowed to explore further
                    print("starting exploring node ", prefix_rpid)
                    reactants_l = REACTION_rec[prefix_rpid][0]["products"] + REACTION_rec[prefix_rpid][0]["unreacted"]
                    print("reactants for this node: ", reactants_l)
                    singlecycle = SingleCycleSearch(STRUC_rec, REACTION_rec, setting_d, wd)
                    explore_dl = singlecycle.explore(reactants_l, prefix_rpid)
                    for reactiontype_d in explore_dl:
                        for count, info_d in reactiontype_d.items():
                            # get products geom
                            info_d = singlecycle.follow_up_treatment_gfn(info_d)
                            # refine products info
                            info_d = singlecycle.follow_up_treatment_energy(info_d)
                            # now evaluate this reaction path
                            reaction_id = prefix_rpid + "_" + str(count)
                            # update rec
                            STRUC_rec, REACTION_rec = singlecycle.evalueate_path(info_d, reaction_id)
                            singlecycle = SingleCycleSearch(STRUC_rec, REACTION_rec, setting_d, wd)
            step += 1
    return STRUC_rec, REACTION_rec

def ReactionReactor(wd, STRUC_rec, setting_d, REACTION_rec, step=0, REFSTRUC_rec={}):
    import random
    setting_inspector = DictionaryInspector(setting_d)
    ncycle = setting_inspector.get_value("searching_setting.general.cycle", 1)
    print("ncycle=", ncycle)
    while step < ncycle:
        print("begin cycle ", step)
        node = [reaction_id for reaction_id in REACTION_rec.keys() if len(reaction_id.split("_")) == step + 1]
        if not node:
            node = ["0"]
        # tmp for now: random pick 2 strucs as reactansts in this cycle
        reactants_l = random.sample(list(STRUC_rec.keys()), 2)
        print("node = ", node)
        for prefix_rpid in node:
            if len(node) ==1:
                status = True
            else:
                status = REACTION_rec[prefix_rpid][-1]
            if status:
                singlecycle = SingleCycleSearch(STRUC_rec, REACTION_rec, setting_d, wd)
                explore_dl = singlecycle.explore(reactants_l, prefix_rpid)
                for reactiontype_d in explore_dl:
                    for count, info_d in reactiontype_d.items():
                        # get products geom
                        info_d = singlecycle.follow_up_treatment_gfn(info_d)
                        # refine products info
                        info_d = singlecycle.follow_up_treatment_energy(info_d)
                        # now evaluate this reaction path
                        reaction_id = prefix_rpid + "_" + str(count)
                        # update rec
                        STRUC_rec, REACTION_rec = singlecycle.evalueate_path(info_d, reaction_id)
                        singlecycle = SingleCycleSearch(STRUC_rec, REACTION_rec, setting_d, wd)
                step += 1
    return STRUC_rec, REACTION_rec


class SingleCycleSearch:
    def __init__(self, STRUC_rec, REACTION_rec, setting_d, wd):
        self.STRUC_rec = STRUC_rec
        self.REACTION_rec = REACTION_rec
        self.setting_d = setting_d
        self.wd = wd

    def ML_prepare(self, struc):
        STRUC_inspector = DictionaryInspector(self.STRUC_rec)
        setting_inspector = DictionaryInspector(self.setting_d)

        # there is a bug in xtb_ml for H
        coord = STRUC_inspector.get_value(f"{struc}.struc_info.general.coord")
        if coord.shape[0] == 1 and coord.iloc[0, 0] == "H":
            wboe = None
            nel = STRUC_inspector.get_value(f"{struc}.struc_info.general.nel")
            ml_df = pd.DataFrame(columns=['index', 'HOAO_a (eV)', 'HOAO_b (eV)', 'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom'])
            if nel == 0:
                # H+
                ml_df.loc[0,:] = [0,-np.inf,-np.inf,-np.inf,-np.inf,1]
                STRUC_inspector.set_value(f"{struc}.struc_info.reactivity.ml", ml_df)
                STRUC_inspector.set_value(f"{struc}.struc_info.reactivity.bo", wboe)
                return
            else:
                # defalut routine
                pass

        # # check if H atom #todo fix xtb ml version for H radical
        # coord = STRUC_inspector.get_value(f"{struc}.struc_info.general.coord")
        # if coord.shape[0]==1 and coord.iloc[0,0] =="H":
        #     wboe = None
        #     nel = STRUC_inspector.get_value(f"{struc}.struc_info.general.nel")
        #     ml_df = pd.DataFrame(columns=['index','HOAO_a (eV)', 'HOAO_b (eV)','LUAO_a (eV)', 'LUAO_b (eV)', 'Atom'])
        #     if nel == 0:
        #         # H+
        #         ml_df.loc[0,:] = [0,-np.inf,-np.inf,-np.inf,-np.inf,1]
        #     elif nel ==1:
        #         ml_df.loc[0, :] = [0, -10.7072, -10.7072, np.inf, np.inf, 1]
        #     else:
        #         ml_df.loc[0, :] = [0, -1.8430, -1.8430, np.inf, np.inf, 1]
        #     STRUC_inspector.set_value(f"{struc}.struc_info.reactivity.ml", ml_df)
        #     STRUC_inspector.set_value(f"{struc}.struc_info.reactivity.bo", wboe)
        #     return

        # do ml
        mlwd = _get_reactive_features(self.wd,struc,self.STRUC_rec,self.setting_d)


        # read ml
        ml_feature = pd.read_csv(os.path.join(mlwd, "ml_feature.csv"), index_col=False)
        ml_feature.insert(0, 'index', range(len(ml_feature)), allow_duplicates=False)
        ml_feature.set_index("index")
        if 'HOAO_a (eV)' in ml_feature.columns:
            ml_df = ml_feature[['index',
                                'HOAO_a (eV)', 'HOAO_b (eV)',
                                'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
        else:
            ml_df = ml_feature[['index',
                                'HOAO (eV)', 'LUAO (eV)', 'Atom']]
            ml_df.rename(columns={'HOAO (eV)': 'HOAO_a (eV)'})
            ml_df.rename(columns={'LUAO (eV)': 'LUAO_a (eV)'})
            ml_df.insert(len(ml_df.columns), "HOAO_b (eV)", ml_df.loc[:, "HOAO_a (eV)"].copy())
            ml_df.insert(len(ml_df.columns), "LUAO_b (eV)", ml_df.loc[:, "LUAO_a (eV)"].copy())

        ###### del H radical beta afo features (no beta electron at all) ######
        nel = STRUC_inspector.get_value(f"{struc}.struc_info.general.nel")
        if nel == 1:
            ml_df.loc[:, "HOAO_b (eV)"] = -np.inf
            ml_df.loc[:, "LUAO_b (eV)"] = np.inf

        STRUC_inspector.set_value(f"{struc}.struc_info.reactivity.ml", ml_df)

        # set bo
        if os.path.exists(os.path.join(mlwd, "wboe")) and os.path.getsize(os.path.join(mlwd, "wboe")) != 0:
            wboe = pd.read_csv(os.path.join(mlwd, "wboe"), header=None, sep='\s+')
        elif os.path.exists(os.path.join(mlwd, "bde")) and os.path.getsize(os.path.join(mlwd, "bde")) != 0:
            wboe = pd.read_csv(os.path.join(mlwd, "bde"), header=None, sep='\s+')
        else:
            wboe = None
        STRUC_inspector.set_value(f"{struc}.struc_info.reactivity.bo", wboe)
        return

    def explore(self, reactants, strnr):
        STRUC_inspector = DictionaryInspector(self.STRUC_rec)
        setting_inspector = DictionaryInspector(self.setting_d)

        # complete ml feature
        for struc in reactants:
            if STRUC_inspector.get_value(f"{struc}.struc_info.reactivity.ml") is None or STRUC_inspector.get_value(
                    f"{struc}.struc_info.reactivity.bo") is None:
                print("prepare ML features for struc ",struc)
                self.ML_prepare(struc)

        # starting searching
        searchcount = 0
        output_interA, output_D, output_intraA = {}, {}, {}
        if setting_inspector.get_value("searching_setting.inter") and len(reactants) > 1:
            ## interA
            merged_setting = {**setting_inspector.get_value("searching_setting.inter"),
                              **setting_inspector.get_value("searching_setting.general")}
            from aRST.script.interA import simulation_interA
            print("starting searching interA posiibilities")
            output_interA = simulation_interA(reactants,
                                              self.STRUC_rec,
                                              merged_setting,
                                              setting_inspector.get_value("xtb_setting"),
                                              self.wd, strnr, count=searchcount)
            searchcount += len(output_interA.keys())

        ## intraD
        if setting_inspector.get_value("searching_setting.disso"):
            merged_setting = {**setting_inspector.get_value("searching_setting.disso"),
                              **setting_inspector.get_value("searching_setting.general")}
            from aRST.script.intraD import simulation_intraD
            print("starting searching Dissociation posiibilities")
            output_D = simulation_intraD(reactants,  # reactants list
                                         self.STRUC_rec,  # correspongding struc info in dict
                                         merged_setting,
                                         setting_inspector.get_value("xtb_setting"),
                                         self.wd, strnr, count=searchcount)
            searchcount += len(output_D.keys())

        # intraA
        if setting_inspector.get_value("searching_setting.intra"):
            from aRST.script.intraA import simulation_interA
            merged_setting = {**setting_inspector.get_value("searching_setting.intra"),
                              **setting_inspector.get_value("searching_setting.general")}
            print("starting searching intraA posiibilities")
            output_intraA = simulation_interA(reactants,  # reactants list
                                              self.STRUC_rec,  # correspongding struc info in dict
                                              merged_setting,
                                              setting_inspector.get_value("xtb_setting"),
                                              self.wd, strnr, count=searchcount)
            searchcount += len(output_intraA.keys())
        explore_l = [output_interA, output_D, output_intraA]
        return explore_l

    def follow_up_treatment_gfn(self, tmp_path_d):
        import random
        import string
        # doing sep opt for inter_fragments
        # output_d[count] = ({reactants info},{endxyz info})
        endxyz_info = tmp_path_d[1]
        fragment_l = seperate_fragment(endxyz_info["coord"], endxyz_info["bo"], endxyz_info["atomic_charge"])
        tmp_path_d[1] = []
        if len(fragment_l) == 1:
            tmp_path_d[1].append({"coord": endxyz_info["coord"],
                                  "atomic_charge": endxyz_info["atomic_charge"],
                                  "charge": endxyz_info["charge"],
                                  "nel": endxyz_info["nel"],
                                  "atomic_connection": endxyz_info["atomic_connection"]})
            return tmp_path_d
        else:
            setting_inspector = DictionaryInspector(self.setting_d)
            # set a random folder to do opt
            for frindex in range(len(fragment_l)):
                coord = fragment_l[frindex]["coord"]
                if coord.shape[0] > 1:
                    random_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
                    optwd = os.path.join(self.wd, "inter_fragments", random_name, "xtb")
                    create_folder(optwd)
                    write_xyz(coord, "coord")

                    if fragment_l[frindex]["nel"] % 2 == 0:
                        uhf = 0
                    else:
                        uhf = 1
                    status = CallxTB(optwd, {"method": setting_inspector.get_value("xtb_setting.method", "gfn2"),
                                             "coordfn": "coord",
                                             "optlevel": setting_inspector.get_value("xtb_setting.optlevel", None),
                                             "charge": fragment_l[frindex]["charge"],
                                             "uhf": uhf,
                                             "alpb": setting_inspector.get_value("xtb_setting.alpb", "h2o"),
                                             "etemp": setting_inspector.get_value("xtb_setting.struc_etemp", "4000")
                                             }).opt()
                    if status == 0:
                        # # check here!!!!
                        # print(fragment_l[frindex])
                        # print(optwd)
                        constrain_l, ele = read_topomol(os.path.join(optwd, 'xtbtopo.mol')) if os.path.exists(
                            os.path.join(optwd, 'xtbtopo.mol')) else [],[]
                        tmp_path_d[1].append({"coord": pd.read_table(os.path.join(optwd, 'xtbopt.xyz'),
                                                                     header=None, skiprows=2, sep='\s+'),
                                              "atomic_charge": pd.read_table(os.path.join(optwd, "charges"),
                                                                             header=None, sep='\s+'),
                                              "charge": fragment_l[frindex]["charge"],
                                              "nel": fragment_l[frindex]["nel"],
                                              "strucwd": random_name,
                                              "atomic_connection": (constrain_l, ele)})
                else:
                    tmp_path_d[1].append({"coord": coord,
                                          "atomic_charge": fragment_l[frindex]["atomic_charge"],
                                          "charge": fragment_l[frindex]["charge"],
                                          "nel": fragment_l[frindex]["nel"],
                                          "atomic_connection": ([], [str(coord.iloc[0, 0])])})
            return tmp_path_d

    def follow_up_treatment_energy(self, tmp_path_d):
        import random
        import string
        info_l = tmp_path_d[1]
        setting_inspector = DictionaryInspector(self.setting_d)
        dftopt = setting_inspector.get_value("orca_setting.dftopt", False)
        functional = setting_inspector.get_value("orca_setting.functional", "r2scan-3c")
        basis = setting_inspector.get_value("orca_setting.basis", "")
        for frindex in range(len(info_l)):
            frinfo_d = info_l[frindex]
            coord = frinfo_d["coord"]
            if coord.shape[0] > 1 and dftopt:
                if frinfo_d["nel"] % 2 == 0:
                    # do both m=1, 3 and choose the sp lowest fr
                    final_m = None
                    sp_dft = np.inf
                    for m in [1, 3]:
                        wd = create_folder(
                            os.path.join(self.wd, "inter_fragments", frinfo_d['strucwd'], "orca", str(m)))
                        write_xyz(coord, "coord")
                        staus = CallORCA(wd, {"functional": functional,
                                              "basis":basis,
                                              "charge": frinfo_d["charge"],
                                              "m": m,
                                              "fn": "coord",
                                              "caltyp": "SP OPT",
                                              "caltyp2": "RKS" if m == 1 else "UKS",
                                              "solmode": setting_inspector.get_value("orca_setting.solmode", ""),
                                              "ncpu": setting_inspector.get_value("orca_setting.ncpu", 1)}).call()
                        tmp_e = get_energy(os.path.join(wd, "ENERGY"))
                        if tmp_e < sp_dft:
                            sp_dft = tmp_e
                            final_m = m
                            tmp_path_d[1][frindex]["coord"] = pd.read_table(os.path.join(wd, 'orca.xyz'),
                                                                            header=None, skiprows=2, sep='\s+')
                            tmp_path_d[1][frindex]["sp_dft"] = tmp_e
                            tmp_path_d[1][frindex]["m"] = final_m
                else:  # do m = 2
                    wd = create_folder(os.path.join(self.wd, "inter_fragments", frinfo_d['strucwd'], "orca", str(2)))
                    write_xyz(coord, "coord")
                    staus = CallORCA(wd, {"functional": functional,
                                          "basis":basis,
                                          "charge": frinfo_d["charge"],
                                          "m": 2,
                                          "fn": "coord",
                                          "caltyp": "SP OPT",
                                          "caltyp2": "UKS",
                                          "solmode": setting_inspector.get_value("orca_setting.solmode", ""),
                                          "ncpu": setting_inspector.get_value("orca_setting.ncpu", 1)}).call()
                    sp_dft = get_energy(os.path.join(wd, "ENERGY"))
                    tmp_path_d[1][frindex]["sp_dft"] = sp_dft
                    tmp_path_d[1][frindex]["m"] = 2

            # complet sp if needed and do trv and solv
            if "sp_dft" not in frinfo_d.keys():
                if frinfo_d["nel"] % 2 == 0:
                    # do both m=1, 3 and choose the sp lowest fr
                    sp_dft = np.inf
                    final_m = None
                    for m in [1, 3]:
                        if "strucwd" in frinfo_d.keys():
                            wd = create_folder(
                                os.path.join(self.wd, "inter_fragments", frinfo_d['strucwd'], "orca", str(m)))
                        else:
                            random_name = ''.join(
                                random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
                            tmp_path_d[1][frindex]["strucwd"] = random_name
                            wd = create_folder(os.path.join(self.wd, "inter_fragments", random_name, "orca", str(m)))

                        write_xyz(coord, "coord")
                        staus = CallORCA(wd, {"functional": functional,
                                              "basis":basis,
                                              "charge": frinfo_d["charge"],
                                              "m": m,
                                              "fn": "coord",
                                              "caltyp": "SP",
                                              "caltyp2": "RKS" if m == 1 else "UKS",
                                              "solmode": setting_inspector.get_value("orca_setting.solmode", ""),
                                              "ncpu": setting_inspector.get_value("orca_setting.ncpu", 1)}).call()
                        tmp_e = get_energy(os.path.join(wd, "ENERGY"))
                        if tmp_e < sp_dft:
                            sp_dft = tmp_e
                            final_m = m
                            tmp_path_d[1][frindex]["sp_dft"] = tmp_e
                            tmp_path_d[1][frindex]["m"] = final_m
                else:
                    if "strucwd" in frinfo_d.keys():
                        wd = create_folder(
                            os.path.join(self.wd, "inter_fragments", frinfo_d['strucwd'], "orca", str(2)))
                    else:
                        random_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
                        tmp_path_d[1][frindex]["strucwd"] = random_name
                        wd = create_folder(os.path.join(self.wd, "inter_fragments", random_name, "orca", str(2)))
                    write_xyz(coord, "coord")
                    staus = CallORCA(wd, {"functional": functional,
                                          "basis":basis,
                                          "charge": frinfo_d["charge"],
                                          "m": 2,
                                          "fn": "coord",
                                          "caltyp": "SP",
                                          "caltyp2": "UKS",
                                          "solmode": setting_inspector.get_value("orca_setting.solmode", ""),
                                          "ncpu": setting_inspector.get_value("orca_setting.ncpu", 1)}).call()
                    sp_dft = get_energy(os.path.join(wd, "ENERGY"))
                    tmp_path_d[1][frindex]["sp_dft"] = sp_dft
                    tmp_path_d[1][frindex]["m"] = 2
            else:
                sp_dft = frinfo_d["sp_dft"]
            # set a tmp value if orca sp error
            if "sp_dft" not in frinfo_d.keys():
                tmp_path_d[1][frindex]["sp_dft"] = np.inf
                tmp_path_d[1][frindex]["m"] = 2 if tmp_path_d[1][frindex]["nel"] % 2 != 0 else 1

            coord = tmp_path_d[1][frindex]["coord"]
            # trv
            if coord.shape[0] > 1:
                wd = create_folder(os.path.join(self.wd, "inter_fragments", frinfo_d['strucwd'], "xtb", "trv"))
                write_xyz(coord, "coord")
                staus = CallxTB(wd, {"method": setting_inspector.get_value("xtb_setting.method", "gfn2"),
                                     "coordfn": "coord",
                                     "optlevel": setting_inspector.get_value("xtb_setting.optlevel", None),
                                     "charge": frinfo_d["charge"],
                                     "uhf": tmp_path_d[1][frindex]["m"] - 1,
                                     "etemp": setting_inspector.get_value("xtb_setting.struc_etemp", "4000"),
                                     "alpb": setting_inspector.get_value("xtb_setting.alpb", "h2o"),
                                     }).trv()
                etrv = get_energy(os.path.join(wd, 'GTRV')) if staus == 0 else 0
            else:
                etrv=0
            tmp_path_d[1][frindex]["trv"] = etrv

            # solv
            wd = create_folder(os.path.join(self.wd, "inter_fragments", frinfo_d['strucwd'], "xtb", "solv"))
            write_xyz(coord, "coord")
            staus = CallxTB(wd, {"method": setting_inspector.get_value("xtb_setting.method", "gfn2"),
                                 "coordfn": "coord",
                                 "optlevel": setting_inspector.get_value("xtb_setting.optlevel", None),
                                 "charge": frinfo_d["charge"],
                                 "uhf": tmp_path_d[1][frindex]["m"] - 1,
                                 "alpb": setting_inspector.get_value("xtb_setting.alpb", "h2o"),
                                 "etemp": setting_inspector.get_value("xtb_setting.struc_etemp", "4000")
                                 }).solv()
            esolv = get_energy(os.path.join(wd, 'GSOLV')) if staus == 0 else 0
            tmp_path_d[1][frindex]["solv"] = esolv
            tmp_path_d[1][frindex]["G_dft"] = esolv + etrv + sp_dft
        return tmp_path_d

    def evalueate_path(self, tmp_path_d, reaction_id):
        mblevel = DictionaryInspector(self.setting_d).get_value("searching_setting.general.mb_level")
        from collections import Counter
        STRUC_inspector = DictionaryInspector(self.STRUC_rec)
        REACTION_inspector = DictionaryInspector(self.REACTION_rec)
        reactants_info = tmp_path_d[0]
        reactants_l = reactants_info["reactants"]
        unreact_l = reactants_info["unreacted"]
        reactants_mbl = [STRUC_inspector.get_value(f"{i}.struc_info.general.mb") for i in reactants_l]
        unreact_mbl = [STRUC_inspector.get_value(f"{i}.struc_info.general.mb") for i in unreact_l]
        products_mbl = [get_molbarcode(tmp_path_d[1][frindex]["coord"], mblevel) for frindex in
                        range(len(tmp_path_d[1]))]
        # record all successful paths (reactants != products) and
        # marked with stuts Ture: allowed to explore further; False: ended in this run.

        # check if reactants = products
        if not "unknow" in products_mbl:
            if Counter(products_mbl) == Counter(reactants_mbl):
                REACTION_inspector.set_value(reaction_id, [{"reactants": reactants_l,
                                                            "unreacted": unreact_l,
                                                            "products": reactants_l,
                                                            "delta_G": 0},
                                                           False])

                return self.STRUC_rec, self.REACTION_rec

        status = True

        # check if reaction exist
        if not "unknow" in products_mbl and not "unknow" in reactants_mbl:
            for refid,item in self.REACTION_rec.items():
                info_d = item[0]
                refreactants = info_d["reactants"]
                refprodcts = info_d["products"]
                refunreact = info_d["unreacted"]
                ref_rmbl = [STRUC_inspector.get_value(f"{i}.struc_info.general.mb") for i in refreactants]
                ref_pmbl = [STRUC_inspector.get_value(f"{i}.struc_info.general.mb") for i in refprodcts]
                ref_umbl = [STRUC_inspector.get_value(f"{i}.struc_info.general.mb") for i in refunreact]

                if Counter(products_mbl) == Counter(ref_pmbl) and Counter(reactants_mbl) == Counter(ref_rmbl):
                    # test unreactants too
                    if Counter(unreact_mbl) == Counter(ref_umbl):
                        status = False

        # check if reaction energy < thres.
        G0 = sum([STRUC_inspector.get_value(f"{i}.struc_info.energy.G_dft") for i in reactants_l])
        G1 = sum([tmp_path_d[1][frindex]["G_dft"] for frindex in range(len(tmp_path_d[1]))])
        delta_G = (G1 - G0) * 627.503
        if delta_G > DictionaryInspector(self.setting_d).get_value("searching_setting.general.elim",300):
            status = False

        # prepare to record reaction info
        # check according to mb if add new struc from products in struc.rec
        products_l = []
        for frindex in range(len(tmp_path_d[1])):
            pmb = products_mbl[frindex]
            addtriger = True
            if pmb != "unknow":
                for struc in self.STRUC_rec.keys():
                    if STRUC_inspector.get_value(f"{struc}.struc_info.general.mb") == pmb:
                        products_l.append(struc)
                        addtriger = False
            if addtriger:
                newstrucid = self.add_newstruc(tmp_path_d[1][frindex])
                products_l.append(newstrucid)

        # add reaction info into reaction rec
        REACTION_inspector.set_value(reaction_id, [{"reactants": reactants_l,
                                                    "unreacted": unreact_l,
                                                    "products": products_l,
                                                    "delta_G": delta_G},
                                                   status])
        return self.STRUC_rec, self.REACTION_rec

    def add_newstruc(self, strucinfo_d):
        inspector_2 = DictionaryInspector(strucinfo_d)  # info of new struc
        inspector_1 = DictionaryInspector(self.STRUC_rec)  # struc database

        newid = str(len(self.STRUC_rec.keys()))

        inspector_1.set_value(f"{newid}.struc_info.general.coord", inspector_2.get_value("coord"))
        inspector_1.set_value(f"{newid}.struc_info.general.charge", inspector_2.get_value("charge"))
        inspector_1.set_value(f"{newid}.struc_info.general.multiplicity", inspector_2.get_value("m"))
        inspector_1.set_value(f"{newid}.struc_info.general.nel", inspector_2.get_value("nel"))
        inspector_1.set_value(f"{newid}.struc_info.general.mb",
                              get_molbarcode(inspector_2.get_value("coord"),
                                             DictionaryInspector(self.setting_d).
                                             get_value("searching_setting.general.mb_level")))
        inspector_1.set_value(f"{newid}.struc_info.general.atomic_connection",
                              inspector_2.get_value("atomic_connection"))
        inspector_1.set_value(f"{newid}.struc_info.energy.sp_dft", inspector_2.get_value("sp_dft"))
        inspector_1.set_value(f"{newid}.struc_info.energy.trv", inspector_2.get_value("trv"))
        inspector_1.set_value(f"{newid}.struc_info.energy.solv", inspector_2.get_value("solv"))
        inspector_1.set_value(f"{newid}.struc_info.energy.G_dft", inspector_2.get_value("G_dft"))

        # move folder
        tmpwd = os.path.join(self.wd, "inter_fragments", inspector_2.get_value("strucwd"))
        tarwd = os.path.join(self.wd, "struc", newid)
        # print(f"move {tmpwd} to {tarwd}")

        if os.path.exists(tmpwd) and not os.path.exists(tarwd):
            shutil.move(tmpwd, tarwd)
        return newid
