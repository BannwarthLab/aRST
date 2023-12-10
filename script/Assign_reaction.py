import os
import shutil
import sys
import traceback

import numpy as np
import pandas as pd
from toolbox import creat_randomfolder, create_folder, write_xyz, get_nelectron, \
    get_charge,get_energy,get_reactionenergy,rm_folderORfile
from head import call_header
from spr import SPR_HL, SPR_dissociation, build_bo_d

from geom import Coord, alignment
from dictionaryInspector import DictionaryInspector
from callsys import CallxTB, CallORCA,CallTerachem,CallsTDA
from struc import StrucOperation
from do_scan_mb import get_endxyz_from_mb


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
    tmpcount:  tmpcount within this reactionpath
    out_struc_d: collect all struc created in this reactionpath
    out_reaction_d: collect all reactions created in this reactionpath

    """
    tmpcount = 1
    out_struc_d = {}
    OUT_Struc = DictionaryInspector(out_struc_d)
    out_reaction_d = {}
    OUT_Reaction = DictionaryInspector(out_reaction_d)

    def __init__(self, wd, STRUC, REACTION, **kwargs):
        self.wd = wd
        self.kwargs = kwargs
        self.struc_inspector = STRUC
        self.reaction_inspector = REACTION


    def do_reaction(self, reactants, reactionpath,current_nstruc, current_nreaction):
        """

        :param current_nreaction:
        :param current_nstruc:
        :param reactants: reaction_idl for this cycle path guess
        :param reactionpath: prefix str of reaction path

        """

        print(f"before guessing for path {reactionpath}: reactants are: {reactants}")




        # results_reactionid_l = []
        searching_setting = self.kwargs.get("searching_setting")
        #
        # prepare reactive features
        self.ML_prepare(reactants)
        react_os, react_ss, HOAO_a, HOAO_b, LUAO_a, LUAO_b = SPR_HL(reactants,
                                                                    self.struc_inspector,
                                                                    allow_hydrogen=searching_setting.get("hydrogen"))
        BO = SPR_dissociation(reactants,
                              self.struc_inspector,
                              allow_hydrogen=searching_setting.get("hydrogen"),type="wboe",
                              allow_halogen=searching_setting.get("halogen"))

        bo = build_bo_d(reactants,
                        self.struc_inspector,
                        allow_hydrogen=searching_setting.get("hydrogen"))
        #
        # start guess
        out_intraa, out_disso, out_intera = [],[],[]



        ############ intramolecuar association ######################
        if bo and (react_os or react_ss):
            out_intraa = self.intra_asso(react_os, react_ss, bo,reactionpath,current_nstruc, current_nreaction)

        ############# dissociation ##################
        if BO:
            out_disso = self.disso(BO,reactionpath,current_nstruc, current_nreaction)


        ##### intermolecular association ############
        if BO and (react_os or react_ss):
            out_intera = self.inter_asso(react_os, react_ss, BO,reactionpath,current_nstruc, current_nreaction)

        """
        out_d:{
            reactionID:[[reactants],[products]]
            }
        out2_d = {reactionID:full_reactants for this cycle}
        """
        # fill in the reactants struc info
        for strucid in reactants:
            if strucid not in self.OUT_Struc.data.keys():
                self.OUT_Struc.data.update({strucid: self.struc_inspector.data[strucid]})

        out_d = {**out_intraa, **out_disso, **out_intera}
        out2_d = {}
        for reaction_id,item in out_d.items():
            out2_d[reaction_id] = reactants
            reactants_l = item[0]
            products_l = item[1]
            unreactants_l = [k for k in reactants if k not in reactants_l]
            str_pathnr = self.OUT_Reaction.get_value(f"{reaction_id}.reaction_path")[0]

            # fill in reaction info
            self.OUT_Reaction.set_reaction_values(self.OUT_Struc, reactants_l,
                                                        products_l, unreactants_l, reaction_id, str_pathnr)
            for strucid in reactants:
                oldvalue = self.OUT_Struc.get_value(f"{strucid}.reaction_path.as_reactant")
                if oldvalue is None:
                    oldvalue = []
                oldvalue.append(str_pathnr)
                self.OUT_Struc.set_value(f"{strucid}.reaction_path.as_reactant",list(set(oldvalue)))
            for strucid in products_l:
                oldvalue = self.OUT_Struc.get_value(f"{strucid}.reaction_path.as_product")
                if oldvalue is None:
                    oldvalue = []
                oldvalue.append(str_pathnr)
                self.OUT_Struc.set_value(f"{strucid}.reaction_path.as_product",list(set(oldvalue)))


        # direct return without del anything

        return self.OUT_Struc.data,self.OUT_Reaction.data,out_d,out2_d

        # # Evaluate reaction routes for repetition of old cycles
        # for out_d in [out_intraa,out_disso,out_intera]:
        #     """
        #     out_d:{
        #         reactionID:[[reactants],[products]]
        #         }
        #     """
        #
        #     for reaction_id,item in out_d.items():
        #         reactants_l = item[0]
        #         products_l = item[1]
        #         unreactants_l = [k for k in reactants if k not in reactants_l]
        #
        #
        #         str_pathnr = reaction_inspector.get_value(f"{reaction_id}.reaction_path")[0]
        #
        #
        #         # test if this reaction repeat
        #         tmpid = reaction_inspector.test_if_reaction_exist()
        #         if tmpid is not None:
        #             # kill temporarily occupied id, but not change others reaction id
        #             del reaction_inspector.data[reaction_id]
        #             old_reactionpath_l = reaction_inspector.get_value(f"{tmpid}.reaction_path")
        #             old_reactionpath_l.append(str_pathnr)
        #             reaction_inspector.set_value(f"{tmpid}.reaction_path",old_reactionpath_l)
        #             results_reactionid_l.append(tmpid)
        #             print(f"reactionID {reaction_id} is repeated as {tmpid}, already been killed")
        #         # else:
        #         #     # reset reactionid
        #         #     newreaction_id = str(len([k for k in reaction_inspector.data.keys()
        #         #                               if reaction_inspector.data[k].get("struc_info", None)
        #         #                               is not None]) + 1).zfill(4)
        #         #     # print("debug reset reactionid")
        #         #     # print(reactants_l, products_l, newreaction_id,str_pathnr)
        #         #     reaction_inspector.set_reaction_values(
        #         #         struc_inspector, reactants_l, products_l,unreactants_l, newreaction_id,str_pathnr)
        #         #     results_reactionid_l.append(newreaction_id)
        #         else:
        #             reaction_inspector.set_reaction_values(struc_inspector, reactants_l, products_l,
        #                                                         unreactants_l, reaction_id,str_pathnr)
        #             results_reactionid_l.append(reaction_id)
        #
        #
        #
        #         # update struc_inspector
        #         for struc_ids, path_attribute in [(reactants_l, "as_reactant"), (products_l, "as_product")]:
        #             for strucid in struc_ids:
        #                 tmp_l = struc_inspector.get_value(f"{strucid}.reaction_path.{path_attribute}")
        #                 tmp_l = [] if tmp_l is None else tmp_l
        #                 tmp_l.append(str_pathnr)
        #                 struc_inspector.set_value(f"{strucid}.reaction_path.{path_attribute}", list(set(tmp_l)))

        # # completing reactivity features to avoid repeating computing in the possibile next cycle
        # try:
        #     self.ML_prepare(list(struc_inspector.data.keys()))
        # except:
        #     print(os.getcwd())
        #     traceback.print_exc()


        # return out_intraa,out_disso,out_intera
    def ML_prepare(self, reactants):

        xtb_setting = self.kwargs.get("xtb_setting")
        struc_wd = os.path.join(self.wd, "struc")
        level = xtb_setting.get('xtblevel', 'gfn2')
        mllevel = xtb_setting.get('mllevel', None)
        alpb = xtb_setting.get('alpb', 'water')
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
                    rc = struc_callxtb.callml(xtbwd,mllevel)
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
                    if rc == 0 :
                        moldenfn = "molden.input"
                        create_folder(os.path.join(struc_wd, f"{strucid}/ml", "stda"))
                        shutil.copyfile(os.path.join(struc_wd, f"{strucid}/ml", moldenfn),
                                        os.path.join(struc_wd, f"{strucid}/ml", "stda", moldenfn))
                        os.chdir(os.path.join(struc_wd, f"{strucid}/ml", "stda"))
                        rc = CallsTDA(moldenfn).callml(stdawd,0)
                        if rc ==0:
                            ml_feature = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml","stda", "ml_feature.csv"))
                            ml_df = ml_feature[['HOAO_a (eV)', 'HOAO_b (eV)',
                                                'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
                            ml_df.reset_index(drop=True, inplace=True)
                            self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.ml", ml_df)
                            if os.path.getsize(os.path.join(struc_wd, f"{strucid}/ml","stda", "bde")) != 0:
                                wboe = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml","stda", "bde"), header=None,
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
                    struc_calltc = CallTerachem("coord",method,charge,m)
                    rc = struc_calltc.get_moldeninput()
                    if rc ==0:
                        moldenfn = "coord.molden"
                        create_folder(os.path.join(struc_wd, f"{strucid}/ml","stda"))
                        shutil.copyfile(os.path.join(struc_wd, f"{strucid}/ml","scr.coord",moldenfn),
                                        os.path.join(struc_wd, f"{strucid}/ml","stda",moldenfn))
                        os.chdir(os.path.join(struc_wd, f"{strucid}/ml", "stda"))
                        rc = CallsTDA(moldenfn).callml(stdawd,3)
                        if rc ==0:
                            ml_feature = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml","stda", "ml_feature.csv"))
                            ml_df = ml_feature[['HOAO_a (eV)', 'HOAO_b (eV)',
                                                'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
                            ml_df.reset_index(drop=True, inplace=True)
                            self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.ml", ml_df)
                            if os.path.getsize(os.path.join(struc_wd, f"{strucid}/ml","stda", "bde")) != 0:
                                wboe = pd.read_csv(os.path.join(struc_wd, f"{strucid}/ml","stda", "bde"), header=None,
                                                   sep='\s+')
                            else:
                                wboe = None
                            self.struc_inspector.set_value(f"{strucid}.struc_info.reactivity.bo", wboe)
                        else:
                            print(f'{struc_wd}/{strucid}/ml wrong')
                    else:
                        print(f'{struc_wd}/{strucid}/ml wrong')

        return

    def intra_asso(self,react_os, react_ss,bo,reactionpath,current_nstruc, current_nreaction):

        searching_setting = self.kwargs.get("searching_setting")
        xtb_setting = self.kwargs.get("xtb_setting")
        level = xtb_setting.get('xtblevel', 'gfn2')
        alpb = xtb_setting.get('alpb', 'water')
        etemp = xtb_setting.get("etemp", str(300))
        scanmode = xtb_setting.get('scanmode', None) # for molbar
        output_d = {}


        # requirement: intramolecular, HL gap low, no bonded, no same atom
        nintraa = searching_setting.get("nintra_association", searching_setting.get("searching_pairs", 1))
        print(f"nintra for reactionpath {reactionpath} : {nintraa}")
        intraa = {}
        for react_d in [react_os, react_ss]:
            for current_key, current_item in react_d.items():
                #  5.93 {('1', '1'): [['0', '0']]}
                for struc_combi, at_l in current_item.items():
                    if len(set(struc_combi)) == 1:  # intramolecular
                        struc = struc_combi[0]
                        if 1 not in [len(set(k)) for k in at_l]:  # not same atom
                            for pair in at_l:
                                if pair in bo[struc] or pair[::-1] in bo[struc]:  # not bonded
                                    continue
                                else:
                                    intraa[current_key] = struc, at_l
        intraa = {key: value for key, value in sorted(intraa.items())[:nintraa]}
        # call scan
        if intraa:
            for gap, items in intraa.items():
                # set tmp reaction id
                reaction_id = str(current_nreaction.value+1).zfill(4)
                current_nreaction.value +=1
                self.OUT_Reaction.set_value(f"{reaction_id}",{}) # occupy a position, could be killed after evaluation

                strucid = intraa[gap][0].zfill(4)
                at_l = intraa[gap][1]
                scan_wd = os.path.join(self.wd, "scan", f"{reaction_id}")
                create_folder(scan_wd)
                coord = self.struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
                write_xyz(coord, "coord")
                charge = self.struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
                uhf = self.struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1
                # check dihedral to rerange at_l
                if len(at_l) == 2:
                    at1, at2, at3, at4 = int(at_l[0][0]), int(at_l[0][1]), int(at_l[1][1]), int(at_l[1][0])
                    if Coord(coord).get_dihedral_angle(at1, at2, at3, at4) > 0:
                        at_l = [[at_l[0][0], at_l[1][1]], [at_l[1][0], at_l[0][1]]]
                push_l = [(str(int(k[0]) + 1),  # at1 (index starts with 1)
                           str(int(k[1]) + 1),  # at2 (index starts with 1)
                           str(np.round(Coord(coord).get_ats_distance(int(k[0]), int(k[1])), 3)),  # d_12
                           str(Coord(coord).get_rr_distance(int(k[0])) +
                               Coord(coord).get_rr_distance(int(k[1])))  # rr_12
                           ) for k in at_l]
                if scanmode is None or scanmode == "xtb":
                    # call xtbscan and record reactionid
                    scan_setting = {
                        "constrain": {
                            "k": xtb_setting.get("constrain_strength", str(1)),
                            "distance": None
                        },
                        "scan": {
                            "mode": "concerted" if len(at_l) == 2 else None,
                            "distance": push_l if push_l else None,
                            "step": str(30)
                        }
                    }
                    struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
                    rc = struc_callxtb.call_scan(scan_setting)
                    if rc == 0:
                        # scan success
                        str_pathnr = reactionpath + str(self.tmpcount)
                        self.tmpcount += 1
                        self.OUT_Reaction.set_value(f"{reaction_id}.reaction_path", [str_pathnr])
                        print(f"resigning reaction tmp_id {reaction_id} to reaction path {str_pathnr}")

                        endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                        chargefile = os.path.join(scan_wd, "charges")
                        nel = get_nelectron(os.path.join(scan_wd, "xtb.out"))
                        products_l = self.followup_processsing_of_endxyz(endxyz, uhf + 1, charge, nel,
                                                                         chargefile,current_nstruc)
                        output_d[reaction_id] = [[strucid],products_l]
                else: # use molbar
                    # do gfn2 sp
                    struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
                    rc = struc_callxtb.callsp()
                    if rc== 0:
                        str_pathnr = reactionpath + str(self.tmpcount)
                        self.tmpcount += 1
                        self.OUT_Reaction.set_value(f"{reaction_id}.reaction_path", [str_pathnr])
                        print(f"resigning reaction tmp_id {reaction_id} to reaction path {str_pathnr}")

                        molfile_path = os.path.join(scan_wd,"xtbtopo.mol")
                        endxyz = get_endxyz_from_mb(molfile_path,push_l)
                        write_xyz(endxyz, "endxyz")
                        chargefile = os.path.join(scan_wd, "charges")
                        nel = get_nelectron(os.path.join(scan_wd, "xtb.out"))
                        products_l = self.followup_processsing_of_endxyz(endxyz, uhf + 1, charge, nel,
                                                                         chargefile, current_nstruc)
                        output_d[reaction_id] = [[strucid], products_l]

        return output_d

    def disso(self,BO,reactionpath,current_nstruc, current_nreaction):

        searching_setting = self.kwargs.get("searching_setting")
        xtb_setting = self.kwargs.get("xtb_setting")
        level = xtb_setting.get('xtblevel', 'gfn2')
        alpb = xtb_setting.get('alpb', 'water')
        etemp = xtb_setting.get("etemp", str(300))
        d = searching_setting.get("ddissociation",2)
        scanmode = xtb_setting.get('scanmode', None)  # for molbar
        output_d = {}

        ndissociation = searching_setting.get("ndissociation", searching_setting.get("searching_pairs", 1))
        print(f"ndissociation : {ndissociation}")

        # tmp for Mechnaophre min_path!!!!
        if self.kwargs.get("terminalat") is not None:
            terminalat_setteing = self.kwargs.get("terminalat")
            strucid = "0001"
            at1 = terminalat_setteing.get("at1",None)
            at2 = terminalat_setteing.get("at2",None)

            reaction_id = str(current_nreaction.value + 1).zfill(4)
            current_nreaction.value += 1

            self.OUT_Reaction.set_value(f"{reaction_id}", {})  # occupy a position
            scan_wd = os.path.join(self.wd, "scan", f"{reaction_id}")
            create_folder(scan_wd)

            coord = self.struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
            write_xyz(coord, "coord")
            d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 3)
            charge = self.struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
            uhf = self.struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1

            push_l = [(str(int(at1) + 1),  # at1
                       str(int(at2) + 1),  # at2
                       str(d12),  # d_12
                       str(d * d12))  # 2 * d_12
                      ]
            scan_setting = {
                "constrain": {
                    "k": xtb_setting.get("constrain_strength", str(1)),
                    "distance": None
                },
                "scan": {
                    "mode": None,
                    "distance": push_l if push_l else None,
                    "step": str(30)
                }
            }
            struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
            rc = struc_callxtb.call_scan(scan_setting)
            if rc == 0:
                str_pathnr = reactionpath + str(self.tmpcount)
                self.tmpcount += 1
                self.OUT_Reaction.set_value(f"{reaction_id}.reaction_path", [str_pathnr])
                print(f"resigning reaction id {reaction_id} to reaction path {str_pathnr}")

                endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                chargefile = os.path.join(scan_wd, "charges")
                nel = get_nelectron(os.path.join(scan_wd, "xtb.out"))
                products_l = self.followup_processsing_of_endxyz(endxyz, uhf + 1, charge, nel,
                                                                 chargefile, current_nstruc)
                output_d[reaction_id] = [[strucid], products_l]
            return output_d

        else:
            for dgap in sorted(BO.keys())[:ndissociation]:
                info_l = BO[dgap][0]
                strucid = info_l[0].zfill(4)
                at1, at2 = info_l[1], info_l[2]

                reaction_id = str(current_nreaction.value+1).zfill(4)
                current_nreaction.value +=1

                self.OUT_Reaction.set_value(f"{reaction_id}", {})  # occupy a position
                scan_wd = os.path.join(self.wd, "scan", f"{reaction_id}")
                create_folder(scan_wd)

                coord = self.struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
                write_xyz(coord, "coord")
                d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 3)
                charge = self.struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
                uhf = self.struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1

                push_l = [(str(int(at1) + 1),  # at1
                           str(int(at2) + 1),  # at2
                           str(d12),  # d_12
                           str(d * d12))  # 2 * d_12
                          ]
                if scanmode is None or scanmode == "xtb":
                    scan_setting = {
                        "constrain": {
                            "k": xtb_setting.get("constrain_strength", str(1)),
                            "distance": None
                        },
                        "scan": {
                            "mode": None,
                            "distance": push_l if push_l else None,
                            "step": str(30)
                        }
                    }
                    struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
                    rc = struc_callxtb.call_scan(scan_setting)
                    if rc == 0:
                        str_pathnr = reactionpath + str(self.tmpcount)
                        self.tmpcount += 1
                        self.OUT_Reaction.set_value(f"{reaction_id}.reaction_path", [str_pathnr])
                        print(f"resigning reaction id {reaction_id} to reaction path {str_pathnr}")

                        endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                        chargefile = os.path.join(scan_wd, "charges")
                        nel = get_nelectron(os.path.join(scan_wd, "xtb.out"))
                        products_l = self.followup_processsing_of_endxyz(endxyz, uhf + 1, charge, nel,
                                                                         chargefile,current_nstruc)
                        output_d[reaction_id] = [[strucid],products_l]
                else:
                    str_pathnr = reactionpath + str(self.tmpcount)
                    self.tmpcount += 1
                    self.OUT_Reaction.set_value(f"{reaction_id}.reaction_path", [str_pathnr])
                    print(f"resigning reaction id {reaction_id} to reaction path {str_pathnr}")

                    struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
                    rc = struc_callxtb.callsp()
                    if rc == 0:
                        molfile_path = os.path.join(scan_wd, "xtbtopo.mol")
                        endxyz = get_endxyz_from_mb(molfile_path, push_l)
                        write_xyz(endxyz, "endxyz")
                        chargefile = os.path.join(scan_wd, "charges")
                        nel = get_nelectron(os.path.join(scan_wd, "xtb.out"))
                        products_l = self.followup_processsing_of_endxyz(endxyz, uhf + 1, charge, nel,
                                                                         chargefile, current_nstruc)
                        output_d[reaction_id] = [[strucid], products_l]

            return output_d

    def inter_asso(self,react_os, react_ss,BO,reactionpath,current_nstruc, current_nreaction):

        searching_setting = self.kwargs.get("searching_setting")
        xtb_setting = self.kwargs.get("xtb_setting")
        level = xtb_setting.get('xtblevel', 'gfn2')
        alpb = xtb_setting.get('alpb', 'water')
        etemp = xtb_setting.get("etemp", str(300))
        scanmode = xtb_setting.get('scanmode', None)  # for molbar
        output_d = {}

        intera = {}
        ninter = searching_setting.get("ninterassociation", searching_setting.get("searching_pairs", 1))
        print(f"ninter : {ninter}")
        for react_d in [react_os, react_ss]:
            for current_key, current_item in react_d.items():
                #  5.93 {('1', '1'): [['0', '0']]}
                for struc_combi, at_l in current_item.items():
                    struc1, struc2 = struc_combi
                    intera[current_key] = struc1, struc2, at_l
        intera = {key: value for key, value in sorted(intera.items())[:ninter]}

        def get_weakest_bo(strucid, BO):
            strucid = strucid.rstrip("0")
            tempBO = {key: [value[0][1], value[0][2]] for key, value in BO.items() if value[0][0] == strucid}
            if not tempBO:
                return None, None, None
            lowest_bo = min(tempBO.keys())
            at_l = tempBO[lowest_bo]
            at1, at2 = at_l[0], at_l[1]
            return lowest_bo, at1, at2

        # call scan
        if intera:
            for gap, items in intera.items():
                reaction_id = str(current_nreaction.value+1).zfill(4)
                current_nreaction.value+=1
                self.OUT_Reaction.set_value(f"{reaction_id}", {})  # occupy a position

                struc1, struc2 = intera[gap][0], intera[gap][1]
                coord1 = self.struc_inspector.get_value(f"{struc1.zfill(4)}.struc_info.general.coord")
                coord2 = self.struc_inspector.get_value(f"{struc2.zfill(4)}.struc_info.general.coord")
                nat_coord1 = coord1.shape[0]
                at_l = intera[gap][2]
                at1, at2 = int(at_l[0][0]), int(at_l[0][1])

                scan_wd = os.path.join(self.wd, "scan", f"{reaction_id}")
                create_folder(scan_wd)
                write_xyz(coord1, f"coord_1")
                write_xyz(coord2, f"coord_2")
                alignedcoord = alignment(coord1, coord2, at1, at2, **xtb_setting)
                write_xyz(alignedcoord, f"coord")

                # check if weakest bond bo < 0.9
                tmpr, tmpat1, tmpat2 = None, None, None
                lowest_bo1, tmp1at1, tmp1at2 = get_weakest_bo(struc1, BO)
                tmpr1 = np.round(Coord(coord1).get_ats_distance(int(tmp1at1), int(tmp1at2)), 3) if tmp1at1 is not None and tmp1at2 is not None else None

                lowest_bo2, tmp2at1, tmp2at2 = get_weakest_bo(struc2, BO)
                if tmp2at1 is not None and tmp2at2 is not None:
                    tmpr2 = np.round(Coord(coord2).get_ats_distance(int(tmp2at1), int(tmp2at2)), 3)
                    tmp2at1, tmp2at2 = str(int(tmp2at1) + nat_coord1), str(int(tmp2at2) + nat_coord1)
                else:tmpr2,tmp2at1, tmp2at2 = None, None, None

                lowest_bo1 = np.inf if lowest_bo1 is None else lowest_bo1
                lowest_bo2 = np.inf if lowest_bo2 is None else lowest_bo2

                if lowest_bo2 < 0.9 and lowest_bo2 < lowest_bo1:
                    tmpr = tmpr2
                    tmpat1, tmpat2 = tmp2at1, tmp2at2
                elif lowest_bo1 < 0.9 and lowest_bo1 < lowest_bo2:
                    tmpr = tmpr1
                    tmpat1, tmpat2 = tmp1at1, tmp1at2

                # check if degerated reaction pairs
                if len(at_l) == 2:
                    at1, at2, at3, at4 = int(at_l[0][0]), int(at_l[0][1]) + nat_coord1, \
                                         int(at_l[1][1]), int(at_l[1][0]) + nat_coord1
                    if Coord(alignedcoord).get_dihedral_angle(at1, at2, at3, at4) > 0:
                        at_l = [[at_l[0][0], str(int(at_l[1][1]) + nat_coord1)],
                                [at_l[1][0], str(int(at_l[0][1]) + nat_coord1)]]
                else:
                    at_l = [[at_l[0][0], str(int(at_l[0][1]) + nat_coord1)]]

                # scan setting
                push_l = [(str(int(k[0]) + 1),  # at1 (index starts with 1)
                           str(int(k[1]) + 1),  # at2 (index starts with 1)
                           str(np.round(Coord(alignedcoord).get_ats_distance(int(k[0]), int(k[1])), 3)),  # d_12
                           str((Coord(alignedcoord).get_rr_distance(int(k[0])) +
                                Coord(alignedcoord).get_rr_distance(int(k[1])) * 0.99))  # rr_12
                           ) for k in at_l]

                charge = sum([self.struc_inspector.get_value(f"{k}.struc_info.general.charge")
                              for k in [struc1.zfill(4), struc2.zfill(4)]])
                uhf = sum([self.struc_inspector.get_value(f"{k}.struc_info.general.multiplicity")
                           for k in [struc1.zfill(4), struc2.zfill(4)]]) - 2

                if scanmode is None or scanmode == "xtb":
                    # call xtbscan and record reactionid
                    scan_setting = {
                        "constrain": {
                            "k": xtb_setting.get("constrain_strength", str(1)),
                            "distance": [str(int(tmpat1) + 1),
                                         str(int(tmpat2) + 1),
                                         str(tmpr * 1.5)] if tmpr is not None else None # tmpr is a loosed bond
                        },
                        "scan": {
                            "mode": "concerted" if len(at_l) == 2 else None,
                            "distance": push_l if push_l else None,
                            "step": str(30)
                        }
                    }
                    # call NEUTRAL scan
                    struc_callxtb = CallxTB(level, "coord", 0, uhf, alpb, etemp)
                    rc = struc_callxtb.call_scan(scan_setting)
                    if rc == 0:
                        str_pathnr = reactionpath + str(self.tmpcount)
                        self.tmpcount += 1
                        self.OUT_Reaction.set_value(f"{reaction_id}.reaction_path", [str_pathnr])
                        print(f"resigning reaction id {reaction_id} to reaction path {str_pathnr}")

                        # first constrained opt with charge
                        tmp_endxyz = pd.read_table(os.path.join(scan_wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                        rwd = creat_randomfolder(scan_wd)
                        write_xyz(tmp_endxyz,"coord")
                        copt_setting = {"distance": push_l if push_l else None}
                        extracmd = xtb_setting.get('extracmd', None)

                        rc = CallxTB(level, "coord", charge, uhf, alpb, etemp).callConstrainedOpt(copt_setting,extracmd)

                        endxyz = pd.read_table(os.path.join(rwd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                        chargefile = os.path.join(rwd, "charges")
                        nel = get_nelectron(os.path.join(rwd, "xtb.out"))
                        products_l = self.followup_processsing_of_endxyz(endxyz, uhf + 1, charge, nel,
                                                                         chargefile,current_nstruc)
                        output_d[reaction_id] = [[struc1.zfill(4),struc2.zfill(4)],products_l]
                else:
                    struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
                    rc = struc_callxtb.callsp()
                    if rc == 0:
                        str_pathnr = reactionpath + str(self.tmpcount)
                        self.tmpcount += 1
                        self.OUT_Reaction.set_value(f"{reaction_id}.reaction_path", [str_pathnr])
                        print(f"resigning reaction id {reaction_id} to reaction path {str_pathnr}")

                        molfile_path = os.path.join(scan_wd, "xtbtopo.mol")
                        endxyz = get_endxyz_from_mb(molfile_path, push_l)
                        write_xyz(endxyz, "endxyz")
                        chargefile = os.path.join(scan_wd, "charges")
                        nel = get_nelectron(os.path.join(scan_wd, "xtb.out"))
                        products_l = self.followup_processsing_of_endxyz(endxyz, uhf + 1, charge, nel,
                                                                         chargefile, current_nstruc)
                        output_d[reaction_id] = [[struc1.zfill(4),struc2.zfill(4)],products_l]

        return output_d

    def followup_processsing_of_endxyz(self, endxyz, m, charge, nel, charge_file,current_nstruc):

        xtb_setting = self.kwargs.get("xtb_setting")
        level = xtb_setting.get('xtblevel', 'gfn2')
        alpb = xtb_setting.get('alpb', 'water')
        extracmd = xtb_setting.get('extracmd', None)
        etemp = xtb_setting.get("etemp", str(300))


        fragment_wdl = StrucOperation(endxyz, m, charge).get_nfragement()

        def trement_of_multi_fragments(fragment_wdl):
            tmp_l = []
            nel_tmp = 0
            for i in fragment_wdl:
                coord_i = pd.read_table(os.path.join(i, "coord.xyz"), header=None, skiprows=2, sep='\s+')
                charge_i = get_charge(coord_i, endxyz, charge_file)

                # check struc if already exist in the struc_inspector
                tmp_mb = StrucOperation(coord_i,1,charge_i).get_molbarcode() # charge and uhf not so importent here
                keys = self.struc_inspector.get_keys_from_value(tmp_mb)
                print(f"check struc if already exist: struc in tmpwd {i}, keys = {keys}")
                keys = [k for k in keys if not "ref" in k] # remove refproduct struc

                if not keys: # new struc
                    struc_id = str(current_nstruc.value+1).zfill(4)
                    current_nstruc.value+=1
                    tmp_l.append(struc_id)

                    self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.coord", coord_i)
                    self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.mb", tmp_mb)
                    self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.formula", tmp_mb.split('|')[2])


                    self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.charge", charge_i)
                    wd = create_folder(os.path.join(self.wd, "struc", f"{struc_id}/xtbopt"))
                    write_xyz(coord_i,"coord")
                    rc = CallxTB(level, "coord", charge_i, 0, alpb, etemp).callsp()
                    nel_i = get_nelectron(os.path.join(wd, "xtb.out"))
                    nel_tmp += nel_i
                    if nel_i % 2 == 0:
                        self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.multiplicity", 1)
                    else:
                        self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.multiplicity", 2)
                else:
                    struc_id = keys[0].split(".")[0]
                    wd = os.path.join(self.wd, "struc", f"{struc_id}/xtbopt")
                    nel_i = get_nelectron(os.path.join(wd, "xtb.out"))
                    nel_tmp += nel_i
                    self.OUT_Struc.data.update({struc_id: self.struc_inspector.data[struc_id]})
                    tmp_l.append(struc_id)

            if nel_tmp != nel:
                print(f"ERROR: sum of tmp fragment {tmp_l} nel {nel_tmp} doesn't equare to "
                      f"endxyz in {charge_file} nel {nel}!!!")
            return tmp_l

        # determine basis features
        if len(fragment_wdl) > 1:
            tmp_l = trement_of_multi_fragments(fragment_wdl)
        else:
            tmp_l = []
            # do gfnopt at endxyz to double-check if struc is stable
            rwd = creat_randomfolder(self.wd)
            write_xyz(endxyz, "coord")
            rc = CallxTB(level, "coord", charge, 0, alpb, etemp).callopt(extracmd)
            if rc == 0:
                coordtmp = pd.read_table(os.path.join(rwd, "xtbopt.xyz"), header=None, skiprows=2, sep='\s+')
                test_nfrgments = StrucOperation(coordtmp, m, charge).get_nfragement()
                if len(test_nfrgments) == 1: # struc stable
                    os.chdir(self.wd)
                    # rm_folderORfile(rwd)
                    tmp_mb = StrucOperation(endxyz, 1, charge).get_molbarcode()
                    keys = self.struc_inspector.get_keys_from_value(tmp_mb)
                    keys = [k for k in keys if not "ref" in k]  # remove refproduct struc

                    if not keys: # new struc
                        struc_id = str(current_nstruc.value+1).zfill(4)
                        current_nstruc.value+=1
                        tmp_l.append(struc_id)
                        self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.coord", endxyz)
                        self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.charge", charge)
                        self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.multiplicity", m)
                        self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.mb", tmp_mb)
                        self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.formula", tmp_mb.split('|')[2])
                    else:
                        struc_id = keys[0].split(".")[0]
                        tmp_l.append(struc_id)
                        self.OUT_Struc.data.update({struc_id: self.struc_inspector.data[struc_id]})

                else:
                    # random_name = rwd.split("/")[-1]
                    # tmp = os.path.abspath(os.path.dirname(rwd))
                    # shutil.move(rwd,os.path.join(tmp,"sep"))
                    # modified_test_nfrgments = [address.replace(random_name, 'sep') for address in test_nfrgments]
                    tmp_l = trement_of_multi_fragments(test_nfrgments)


        # do gfn unconstrained opt and record energy
        level2 = self.kwargs.get("orca_setting",None).get("dft_level",None)
        pcm = self.kwargs.get("orca_setting",None).get("pcm",None)
        for struc_id in tmp_l:
            tmp_d = self.OUT_Struc.get_value(f"{struc_id}.struc_info.general")
            coord = tmp_d["coord"]
            charge = tmp_d["charge"]
            uhf = tmp_d["multiplicity"] - 1
            m = tmp_d["multiplicity"]

            # do gfn unconstrained opt
            wd = create_folder(os.path.join(self.wd, "struc", f"{struc_id}/xtbopt"))
            write_xyz(coord, "coord")
            if coord.shape[0] > 1:
                rc = CallxTB(level, "coord", charge, uhf, alpb, etemp).callopt(extracmd)
                if rc == 0:
                    energy = get_energy("ENERGY")
                    self.OUT_Struc.set_value(f"{struc_id}.struc_info.energy.sp_gfn2", energy)
                    xyz = pd.read_table(os.path.join(wd, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
                    self.OUT_Struc.set_value(f"{struc_id}.struc_info.general.coord", xyz)

            # record dft energy
            if not "gfn" in level2:
                wd = create_folder(os.path.join(self.wd, "struc", f"{struc_id}/orca"))
                coord = self.OUT_Struc.get_value(f"{struc_id}.struc_info.general.coord")
                write_xyz(coord, "coord")
                rc = CallORCA(level2, "coord", charge, m, pcm).callcalculation()
                energy = get_energy("ENERGY") if rc == 0 else None
                self.OUT_Struc.set_value(f"{struc_id}.struc_info.energy.sp_dft", energy)

        return tmp_l



def complete_struc_dftopt(wd,info_d,optlevel,pcm):
    strucid = list(info_d.keys())[0]
    tmp_D = DictionaryInspector(info_d)
    coord = tmp_D.get_value(f"{strucid}.struc_info.general.coord")
    charge = tmp_D.get_value(f"{strucid}.struc_info.general.charge")
    m = tmp_D.get_value(f"{strucid}.struc_info.general.multiplicity")
    energy = tmp_D.get_value(f"{strucid}.struc_info.energy.sp_dft") or np.inf

    if coord.shape[0] > 1:
        wd = os.path.join(wd,"struc",strucid,"orca")

        orcawd = create_folder(wd)

        out_m = 1
        out_coord = None
        if m == 1:
            # do m =1 and m = 3 opt
            for current_m in [1,3]:
                rwd = creat_randomfolder(orcawd)
                write_xyz(coord,"coord")
                rc = CallORCA(optlevel, "coord", charge, current_m, pcm).callcalculation(caltyp='sp OPT')
                if rc == 0:
                    current_energy = get_energy("ENERGY") if rc == 0 else None
                    if current_energy < energy:
                        energy = current_energy
                        out_m = current_m
                        out_coord = pd.read_table(os.path.join(rwd, 'orca.xyz'), header=None, skiprows=2, sep='\s+')
        else:
            # open shell
            rwd = creat_randomfolder(orcawd)
            write_xyz(coord, "coord")
            rc = CallORCA(optlevel, "coord", charge, 2, pcm).callcalculation(caltyp='sp OPT')
            if rc == 0:
                energy = get_energy("ENERGY") if rc == 0 else None
                out_m = 2
                out_coord = pd.read_table(os.path.join(rwd, 'orca.xyz'), header=None, skiprows=2, sep='\s+')
    else:
        out_coord = coord
        out_m = m
    return strucid,out_coord,out_m,energy

def complete_trvandsolv(wd,info_d,level,alpb):
    strucid = list(info_d.keys())[0]
    tmp_D = DictionaryInspector(info_d)
    coord = tmp_D.get_value(f"{strucid}.struc_info.general.coord")
    charge = tmp_D.get_value(f"{strucid}.struc_info.general.charge")
    uhf = tmp_D.get_value(f"{strucid}.struc_info.general.multiplicity")-1
    energy_trv,energy_solv = None,None


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

    return strucid,energy_trv,energy_solv


