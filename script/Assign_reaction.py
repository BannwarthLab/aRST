import os
import shutil
import sys


import numpy as np
import pandas as pd
from toolbox import creat_randomfolder, create_folder, write_xyz, get_nelectron, \
    get_charge, get_energy, rm_folderORfile,get_total_charge

from dictionaryInspector import DictionaryInspector
from callsys import CallxTB, CallORCA, CallTerachem, CallsTDA
from struc import StrucOperation



# def do_assign_reactor(react_l, rec, reaction_rec=None, fail_rec=None):
#     if reaction_rec is None:
#         reaction_rec = {}
#     if fail_rec is None:
#         fail_rec = {}
#     reactants_tmp = []  # reactants for this cycle
#     if len(react_l) > 0:
#         while not reactants_tmp:  # assume that +2/-2 charged molecules are not stable and can only dissociate back
#             reactants_tmp = random.sample(react_l, 1) + random.sample(react_l, 1)
#             reactants_mb = [rec[k]['mb'] for k in reactants_tmp]
#             tot_c = float('inf')
#             try:
#                 tot_c = sum([rec[k]['charge'] for k in reactants_tmp])
#             except:
#                 pass
#             # test if reaction already exist #todo
#             if existreaction(fail_rec, reactants_mb, 'r'):
#                 reactants_tmp = []
#             else:
#                 if existreaction(reaction_rec, reactants_mb, 'r'):
#                     reactants_tmp = []
#                 else:
#                     # test if unusual reactants
#                     if tot_c > 1 or tot_c < -1:
#                         fail_rec[len(fail_rec.keys()) + 1] = reactants_mb
#                         reactants_tmp = []
#     else:
#         reactants_tmp = []
#     return reactants_tmp, fail_rec

class Singlestep:
    """
    tmpcount: tmpcount for reaction number within this reactionpath
    out_struc_d: collect all struc created in this reactionpath
    out_reaction_d: collect all reactions created in this reactionpath

    """
    tmpcount = 1
    out_struc_d = {}
    out_reaction_d = {}

    OUT_Struc = DictionaryInspector(out_struc_d) # reaction info in this step
    OUT_Reaction = DictionaryInspector(out_reaction_d) # reaction info in this step

    def __init__(self, wd, STRUC, REACTION, **kwargs):
        self.wd = wd
        self.kwargs = kwargs
        self.struc_inspector = STRUC
        self.reaction_inspector = REACTION

    def do_reaction(self, reactants, reactionpath, current_nstruc, current_nreaction):
        """
        reactants: reactants id list for this step
        reactionpath: rp label (prefix str of reaction path)
        current_nstruc: current number of collected struc
        current_nreaction: current number of collected reaction

        """

        print(f"before guessing for path {reactionpath}: reactants are: {reactants}")

        # Determine sreach range
        searching_setting = self.kwargs.get("searching_setting")
        xtb_setting = self.kwargs.get("xtb_setting")
        sr_pair = searching_setting.get("searching_pairs", None)
        if sr_pair is not None:
            intra_pair,dis_pair,inter_pair = sr_pair,sr_pair,sr_pair
        else:
            intra_pair = searching_setting.get("nintra_association", None)
            dis_pair = searching_setting.get("ndissociation", None)
            inter_pair = searching_setting.get("ninterassociation", None)

        ## prepare reactive features
        self.ML_prepare(reactants) # run ml

        # get atomic connection info
        from toolbox import read_topomol
        for strucid in reactants:
            if self.struc_inspector.get_value(f"{strucid}.struc_info.atomic_connection") is None:
                coord_i = self.struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
                if coord_i.shape[0]>1:
                    molfile = os.path.join(self.wd, "struc", strucid, "ml", "xtbtopo.mol")
                    constrain_l, ele = read_topomol(molfile)
                else:
                    constrain_l=[]
                    ele = [str(coord_i.loc[0,0])]
                self.struc_inspector.set_value(f"{strucid}.struc_info.atomic_connection", (constrain_l, ele))


        out_d = {} # out_d:{reactionID:[[reactants_l],[products_l]]}
        out2_d = {} # out2_d = {reactionID:full_reactants for this cycle}

        ## start guess
        tmpout_0,tmpout_1,tmpout_2={},{},{}
        # Intramolecular Association
        if intra_pair is not None and intra_pair:
            from intraA import simulation_intraA
            tmpout_0 = simulation_intraA(reactants,reactionpath,self.struc_inspector,searching_setting,
                                   xtb_setting,intra_pair,self.wd)
            #tmpout[count] = ([strucid], at_l, gap, endxyz, chargefile, nel)

        # Dissociation
        if dis_pair is not None and dis_pair:
            from intraD import simulation_intraD
            tmpout_1 = simulation_intraD(reactants, reactionpath, self.struc_inspector, searching_setting,
                                       xtb_setting, dis_pair, self.wd)

        # Intermolecular Assosication
        if inter_pair is not None and inter_pair:
            if self.kwargs.get("forcedconstrain_setting") is not None:
                from interA import force_interA
                tmpout_2 = force_interA(reactants, reactionpath, self.struc_inspector,
                                        self.kwargs.get("forcedconstrain_setting"),
                                        xtb_setting, intra_pair, self.wd)
            else:
                from interA import simulation_interA
                tmpout_2 = simulation_interA(reactants, reactionpath, self.struc_inspector, searching_setting,
                                           xtb_setting, inter_pair, self.wd)


        tmpout = {**tmpout_0, **tmpout_1, **tmpout_2}

        fragment_wd = creat_randomfolder(os.path.join(self.wd,"fragment"))

        for i in tmpout.keys():
            reactants_l,at_l, gap, endxyz, chargefile, nel = tmpout[i]
            charge = get_total_charge(chargefile)

            # processing endxyz to get products_l
            products_l = self.followup_processsing_of_endxyz(endxyz,charge,nel,chargefile,current_nstruc,fragment_wd,reactionpath)
            unreactants_l = [k for k in reactants if k not in reactants_l]

            reactionid = str(current_nreaction.value+1).zfill(4)
            current_nreaction.value +=1
            print(f"in {reactionpath},current_nreaction+=1 now is {current_nreaction.value}")
            str_rp = reactionpath+str(self.tmpcount)
            self.tmpcount +=1

            # fill in struc info (isolated variante to golbal variant)
            for id in self.struc_inspector.data.keys():
                if id not in self.OUT_Struc.data.keys():
                    self.OUT_Struc.data.update({id: self.struc_inspector.data[id]})

            # fill in reaction info (isolated variante to golbal variant)
            self.OUT_Reaction.set_reaction_values(self.OUT_Struc, reactants_l,
                                                  products_l, unreactants_l, reactionid, str_rp,at_l,gap)

            # some not so important info
            for strucid in reactants:
                oldvalue = self.OUT_Struc.get_value(f"{strucid}.reaction_path.as_reactant")
                if oldvalue is None:
                    oldvalue = []
                oldvalue.append(str_rp)
                self.OUT_Struc.set_value(f"{strucid}.reaction_path.as_reactant", list(set(oldvalue)))
            for strucid in products_l:
                oldvalue = self.OUT_Struc.get_value(f"{strucid}.reaction_path.as_product")
                if oldvalue is None:
                    oldvalue = []
                oldvalue.append(str_rp)
                self.OUT_Struc.set_value(f"{strucid}.reaction_path.as_product", list(set(oldvalue)))

            out_d[reactionid] = [reactants_l,products_l]
            out2_d[reactionid] = reactants

        rm_folderORfile(fragment_wd)
        return self.OUT_Struc.data, self.OUT_Reaction.data, out_d, out2_d


    def ML_prepare(self, reactants):

        xtb_setting = self.kwargs.get("xtb_setting")
        struc_wd = os.path.join(self.wd, "struc")
        level = xtb_setting.get('xtblevel', 'gfn2')
        mllevel = xtb_setting.get('mllevel', None)
        alpb = xtb_setting.get('alpb', None)
        etemp = xtb_setting.get("etemp", str(300))
        xtbwd = xtb_setting.get('xtb_wd', None)

        # get atomistic features and bo
        for strucid in reactants:
            if self.struc_inspector.get_value(f"{strucid}.struc_info.reactivity.ml") is None:
                coord = self.struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
                charge = self.struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
                m = self.struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity")

                output_folder = os.path.join(struc_wd, f"{strucid}/ml")
                create_folder(output_folder)
                write_xyz(coord, "coord")

                if mllevel is None or "gfn2" in mllevel or "gfn 2" in mllevel:
                    # call xtb

                    struc_callxtb = CallxTB(level, "coord", charge, m - 1, alpb, etemp)
                    rc = struc_callxtb.callml(xtbwd, mllevel)
                    if rc == 0:
                        # set ml_features
                        ml_feature = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml", "ml_feature.csv"))
                        if 'HOAO_a (eV)' in ml_feature.columns:
                            ml_df = ml_feature[['atomic partial charges', 'coordination number',
                                                'HOAO_a (eV)', 'HOAO_b (eV)',
                                                'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
                        else:
                            ml_df = ml_feature[['atomic partial charges', 'coordination number',
                                                'HOAO (eV)', 'LUAO (eV)', 'Atom']]
                        ml_df.reset_index(drop=True, inplace=True)
                        self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.ml", ml_df)
                        # set bo
                        if os.path.getsize(os.path.join(struc_wd, f"{strucid}/ml", "wboe")) != 0:
                            wboe = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml", "wboe"), header=None, sep='\s+')
                        else:
                            wboe = None
                        self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.bo", wboe)
                    else:
                        sys.exit(f'{struc_wd}/{strucid}/ml wrong')
                elif "gfn" in mllevel:
                    # call xtb but later stda
                    stdawd = xtb_setting.get('stda_wd', None)
                    rc = CallxTB(level, "coord", charge, m - 1, alpb, etemp).get_moldeninput()
                    if rc == 0:
                        moldenfn = "molden.input"
                        create_folder(os.path.join(struc_wd, f"{strucid}/ml", "stda"))
                        shutil.copyfile(os.path.join(struc_wd, f"{strucid}/ml", moldenfn),
                                        os.path.join(struc_wd, f"{strucid}/ml", "stda", moldenfn))
                        os.chdir(os.path.join(struc_wd, f"{strucid}/ml", "stda"))
                        rc = CallsTDA(moldenfn).callml(stdawd, 0)
                        if rc == 0:
                            ml_feature = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml", "stda", "ml_feature.csv"))
                            ml_df = ml_feature[['HOAO_a (eV)', 'HOAO_b (eV)',
                                                'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
                            ml_df.reset_index(drop=True, inplace=True)
                            self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.ml", ml_df)
                            if os.path.getsize(os.path.join(struc_wd, f"{strucid}/ml", "stda", "bde")) != 0:
                                wboe = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml", "stda", "bde"), header=None,
                                                   sep='\s+')
                            else:
                                wboe = None
                            self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.bo", wboe)
                        else:
                            print(f'{struc_wd}/{strucid}/ml wrong')
                    else:
                        print(f'{struc_wd}/{strucid}/ml wrong')

                else:
                    # call terachem and then stda
                    stdawd = xtb_setting.get('stda_wd', None)
                    method = mllevel.split("_")[1]
                    struc_calltc = CallTerachem("coord", method, charge, m)
                    rc = struc_calltc.get_moldeninput()
                    if rc == 0:
                        moldenfn = "coord.molden"
                        create_folder(os.path.join(struc_wd, f"{strucid}/ml", "stda"))
                        shutil.copyfile(os.path.join(struc_wd, f"{strucid}/ml", "scr.coord", moldenfn),
                                        os.path.join(struc_wd, f"{strucid}/ml", "stda", moldenfn))
                        os.chdir(os.path.join(struc_wd, f"{strucid}/ml", "stda"))
                        rc = CallsTDA(moldenfn).callml(stdawd, 3)
                        if rc == 0:
                            ml_feature = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml", "stda", "ml_feature.csv"))
                            ml_df = ml_feature[['HOAO_a (eV)', 'HOAO_b (eV)',
                                                'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
                            ml_df.reset_index(drop=True, inplace=True)
                            self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.ml", ml_df)
                            if os.path.getsize(os.path.join(struc_wd, f"{strucid}/ml", "stda", "bde")) != 0:
                                wboe = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml", "stda", "bde"), header=None,
                                                   sep='\s+')
                            else:
                                wboe = None
                            self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.bo", wboe)
                        else:
                            print(f'{struc_wd}/{strucid}/ml wrong')
                    else:
                        print(f'{struc_wd}/{strucid}/ml wrong')

        return


    def followup_processsing_of_endxyz(self, endxyz,charge, nel, charge_file, current_nstruc,fragment_wd,reactionpath):

        xtb_setting = self.kwargs.get("xtb_setting")
        level = xtb_setting.get('xtblevel', 'gfn2')
        alpb = xtb_setting.get('alpb', None)
        extracmd = xtb_setting.get('extracmd', None)
        etemp = xtb_setting.get("etemp", str(300))

        tmp_products=[]


        def creat_new_struc(coord_i, charge_i, nel_i, mb_i,energy_i):
            struc_id = str(current_nstruc.value + 1).zfill(4)
            current_nstruc.value += 1
            print(f"in {reactionpath},current_nstruc+=1 now is {current_nstruc.value}")


            self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.coord", coord_i)
            self.OUT_Struc.set_value(f"{struc_id}.struc_info.energy.sp_gfn2", energy_i)
            self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.charge", charge_i)
            self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.nel", nel_i)
            self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.mb", mb_i)
            if "|" in mb_i:
                self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.formula", mb_i.split('|')[2])
            else:
                self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.formula", mb_i)

            if nel_i%2==1:
                self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.multiplicity", 2)
            else:self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.multiplicity", 1)
            return struc_id

        def trement_of_multi_fragments(fragment_wdl):
            nel_tmp = []
            tmp_STRUC={}
            for i in fragment_wdl:
                tmp_STRUC[i]={}
                coord_i = pd.read_table(os.path.join(i, "coord.xyz"), header=None, skiprows=2, sep='\s+')
                charge_i = get_charge(coord_i, endxyz, charge_file)
                tmp_STRUC[i]["coord"] = coord_i
                tmp_STRUC[i]["c"] = charge_i

                # do sp or unconstrained gfnopt for fragment

                rwd = creat_randomfolder(i)
                write_xyz(coord_i, "coord")
                if coord_i.shape[0] > 1:
                    rc = CallxTB(level, "coord", charge_i, None, alpb, etemp).callopt(extracmd)
                    coord_i = pd.read_table(os.path.join(rwd, "xtbopt.xyz"), header=None, skiprows=2, sep='\s+')
                    energy_i = get_energy(os.path.join(rwd, "ENERGY"))
                    tmp_STRUC[i]["coord"] = coord_i
                    tmp_STRUC[i]["energy"] = energy_i


                else:
                    rc = CallxTB(level, "coord", charge_i, None, alpb, etemp).callsp()
                    energy_i = get_energy(os.path.join(rwd, "ENERGY"))
                    tmp_STRUC[i]["energy"] = energy_i
                nel_i = get_nelectron(os.path.join(rwd, "xtb.out"))
                tmp_STRUC[i]["nel"] = nel_i
                nel_tmp.append(nel_i)
                rm_folderORfile(rwd)

            if sum(nel_tmp) != nel:
                if abs(sum(nel_tmp)-nel) < 2:
                    # adujest the charge in bigest struc
                    maxindex = nel_tmp.index(max(nel_tmp))
                    maxkey = list(tmp_STRUC.keys())[maxindex]
                    tmp_STRUC[maxkey]["c"] = tmp_STRUC[maxkey]["c"] + \
                                                                  sum(nel_tmp)-nel

            # now record fragments
            for i in fragment_wdl:
                coord_i,charge_i,nel_i,energy_i = tmp_STRUC[i]["coord"],tmp_STRUC[i]["c"],\
                                                  tmp_STRUC[i]["nel"],tmp_STRUC[i]["energy"]
                mb_i = StrucOperation(coord_i, None, None).get_molbarcode(opt=False)
                if mb_i is not None:
                    keys = self.struc_inspector.find_strucid_from_mb(mb_i)
                    keys = [k for k in keys if not "ref" in k]  # remove refproduct struc
                else:
                    print(f"Morbar is error with {coord_i}, treate it as new struc")
                    keys = []

                if not keys: # newstruc
                    struc_id = creat_new_struc(coord_i,charge_i,nel_i,mb_i,energy_i)
                    tmp_products.append(struc_id)
                else:
                    trigger2 = True
                    tmp = []


                    for id in keys:
                        tmpcharge = self.struc_inspector.get_value(f"{id}.struc_info.general.charge")
                        tmpnel = self.struc_inspector.get_value(f"{id}.struc_info.general.nel")
                        if tmpcharge == charge_i and tmpnel == nel_i:
                            trigger2 = False
                            # exist strucid
                            tmp.append(id)
                    if trigger2:
                        struc_id = creat_new_struc(coord_i, charge_i, nel_i, mb_i,energy_i)
                        tmp_products.append(struc_id)
                    else:
                        tmp_products.append(tmp[0])  # falls same struc with different id
                        self.OUT_Struc.data.update({tmp[0]: self.struc_inspector.data[tmp[0]]})

            return tmp_products


        fragment_wdl = StrucOperation(endxyz, 0, charge).get_nfragement(fragment_wd)



        # determine basis features
        if len(fragment_wdl) > 1:
            tmp_products = trement_of_multi_fragments(fragment_wdl)
        else:

            # do gfnopt at endxyz to double-check if struc is stable
            trigger = False
            rwd = creat_randomfolder(fragment_wd)
            write_xyz(endxyz, "coord")

            if nel%2==1:
                rc = CallxTB(level, "coord", charge, 1, alpb, etemp).callopt(extracmd)
                coordtmp = pd.read_table(os.path.join(rwd, "xtbopt.xyz"), header=None, skiprows=2, sep='\s+')
                energy = get_energy(os.path.join(rwd, "ENERGY"))
                test_nfrgments = StrucOperation(coordtmp, None, None).get_nfragement(fragment_wd)
                if len(test_nfrgments) == 1:
                    trigger = True
            else:
                for u in [0,2]:
                    CallxTB(level, "coord", charge, u, alpb, etemp).callopt(extracmd)
                    coordtmp = pd.read_table(os.path.join(rwd, "xtbopt.xyz"), header=None, skiprows=2, sep='\s+')
                    energy = get_energy(os.path.join(rwd, "ENERGY"))
                    test_nfrgments = StrucOperation(coordtmp, 1, charge).get_nfragement(fragment_wd)
                    if len(test_nfrgments) == 1:
                        trigger = True
                        break
            if trigger:
                # struc stable
                rm_folderORfile(rwd)
                tmp_mb = StrucOperation(coordtmp, None, None).get_molbarcode(opt=False)
                keys = self.struc_inspector.find_strucid_from_mb(tmp_mb)
                keys = [k for k in keys if not "ref" in k]  # remove refproduct struc

                if not keys:  # new struc
                    struc_id = creat_new_struc(coordtmp, charge, nel,tmp_mb,energy)
                    tmp_products.append(struc_id)
                else:
                    trigger2=True
                    tmp=[]
                    for id in keys:
                        tmpcharge = self.struc_inspector.get_value(f"{id}.struc_info.general.charge")
                        tmpnel = self.struc_inspector.get_value(f"{id}.struc_info.general.nel")
                        if tmpcharge == charge and tmpnel == nel:
                            trigger2 = False
                            # exist strucid
                            tmp.append(id)
                    tmp_products.append(tmp[0]) # falls same struc with different id
                    self.OUT_Struc.data.update({tmp[0]: self.struc_inspector.data[tmp[0]]})

                    if trigger2:
                        struc_id = creat_new_struc(coordtmp, charge,nel, tmp_mb,energy)
                        tmp_products.append(struc_id)
            else: tmp_products = trement_of_multi_fragments(test_nfrgments)

        ######energy part#######

        # do DFT sp and record energy
        level2 = self.kwargs.get("orca_setting", None).get("dft_level", None)
        pcm = self.kwargs.get("orca_setting", None).get("pcm", None)

        for struc_id in tmp_products:
            tmp_d = self.OUT_Struc.get_value(f"{struc_id}.struc_info.general")
            coord = tmp_d["coord"]
            charge = tmp_d["charge"]
            nel = tmp_d["nel"]


            # record dft energy
            if not "gfn" in level2:
                wd = create_folder(os.path.join(self.wd, "struc", f"{struc_id}/orca"))

                write_xyz(coord, "coord")
                if nel%2 ==1:
                    rc = CallORCA(level2, "coord", charge, 2, pcm).callcalculation()
                    self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.m", 2)
                else:
                    rc = CallORCA(level2, "coord", charge, 1, pcm).callcalculation()
                    self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.m", 0)
                energy = get_energy("ENERGY") if rc == 0 else None
                self.OUT_Struc.set_value(f"{struc_id}.struc_info.energy.sp_dft", energy)

        return tmp_products


def complete_struc_dftopt(wd, info_d, optlevel, pcm):
    strucid = list(info_d.keys())[0]
    tmp_D = DictionaryInspector(info_d)
    coord = tmp_D.get_value(f"{strucid}.struc_info.general.coord")
    charge = tmp_D.get_value(f"{strucid}.struc_info.general.charge")
    m = tmp_D.get_value(f"{strucid}.struc_info.general.multiplicity")
    energy = tmp_D.get_value(f"{strucid}.struc_info.energy.sp_dft") or np.inf

    wd = os.path.join(wd, "struc", strucid, "orca")

    orcawd = create_folder(wd)

    out_m = 1
    out_coord = coord

    if m == 1:
        # do m =1 and m = 3 opt
        for current_m in [1, 3]:
            rwd = creat_randomfolder(orcawd)
            write_xyz(coord, "coord")
            rc = CallORCA(optlevel, "coord", charge, current_m, pcm). \
                callcalculation(caltyp='sp OPT' if coord.shape[0] > 1 else 'sp')
            if rc == 0:
                current_energy = get_energy("ENERGY") if rc == 0 else None
                if current_energy < energy:
                    energy = current_energy
                    out_m = current_m
                    out_coord = pd.read_table(os.path.join(rwd, 'orca.xyz'),
                                              header=None, skiprows=2, sep='\s+') if coord.shape[0] > 1 else coord
    else:
        # open shell
        rwd = creat_randomfolder(orcawd)
        write_xyz(coord, "coord")
        rc = CallORCA(optlevel, "coord", charge, 2, pcm). \
            callcalculation(caltyp='sp OPT' if coord.shape[0] > 1 else 'sp')
        if rc == 0:
            energy = get_energy("ENERGY") if rc == 0 else None
            out_m = 2
            out_coord = pd.read_table(os.path.join(rwd, 'orca.xyz'),
                                      header=None, skiprows=2, sep='\s+') if coord.shape[0] > 1 else coord
    return strucid, out_coord, out_m, energy


def complete_trvandsolv(wd, info_d, level, alpb):
    strucid = list(info_d.keys())[0]
    tmp_D = DictionaryInspector(info_d)
    coord = tmp_D.get_value(f"{strucid}.struc_info.general.coord")
    charge = tmp_D.get_value(f"{strucid}.struc_info.general.charge")
    uhf = tmp_D.get_value(f"{strucid}.struc_info.general.multiplicity") - 1
    energy_trv, energy_solv = None, None

    wd1 = os.path.join(wd, "struc", strucid, "trv")
    create_folder(wd1)
    write_xyz(coord, "coord")
    rc = CallxTB(level, "coord", charge, uhf, alpb, None).calltrv()
    if rc == 0:
        energy_trv = get_energy("GTRV")

    wd2 = os.path.join(wd, "struc", strucid, "solv")
    create_folder(wd2)
    write_xyz(coord, "coord")
    rc = CallxTB(level, "coord", charge, uhf, alpb, None).callsolv()
    if rc == 0:
        energy_solv = get_energy("GSOLV")

    return strucid, energy_trv, energy_solv
