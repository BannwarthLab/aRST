import os
import sys
import time
import pandas as pd
from callsys import call_orca_dft
from Subprocess import subprocess, check_xtb_scan, elvauate_reaction
from read import update_dict
import random
from PushAndPull import push, pull, association_pairs, dissociation_pair
from struc import Struc
from Subprocess import existreaction
from operation import collect_sp_xyz, collect_g_xyz
from head import call_header


def do_assign_reactor(react_l, rec, reaction_rec=None, fail_rec=None):
    if reaction_rec is None:
        reaction_rec = {}
    if fail_rec is None:
        fail_rec = {}
    reactants_tmp = []  # reactants for this cycle
    if len(react_l) > 0:
        while not reactants_tmp:  # assume that +2/-2 charged molecules are not stable and can only dissociate back
            reactants_tmp = random.sample(react_l, 1) + random.sample(react_l, 1)
            reactants_mb = [rec[k]['mb'] for k in reactants_tmp]
            tot_c = float('inf')
            try:
                tot_c = sum([rec[k]['charge'] for k in reactants_tmp])
            except:
                pass
            # test if reaction already exist #todo
            if existreaction(fail_rec, reactants_mb, 'r'):
                reactants_tmp = []
            else:
                if existreaction(reaction_rec, reactants_mb, 'r'):
                    reactants_tmp = []
                else:
                    # test if unusual reactants
                    if tot_c > 1 or tot_c < -1:
                        fail_rec[len(fail_rec.keys()) + 1] = reactants_mb
                        reactants_tmp = []
    else:
        reactants_tmp = []
    return reactants_tmp, fail_rec


class Singlestep:
    def __init__(self, wd, **kwargs):
        self.wd = wd
        self.kwargs = kwargs
        self.searching_setting = kwargs.get("searching_setting")
        self.orca_setting = kwargs.get("orca_setting")
        self.xtb_setting = kwargs.get("xtb_setting")

    def do_reaction(self, reactants, stepnr, rec, reaction_rec=None):
        # get atomistic features
        for name in reactants:
            if "ml_feature" not in rec[name].keys():
                coord = rec[name]['coord']
                charg = int(rec[name]['charge'])
                m = int(rec[name]['m'])
                Struc(coord, m, charg).runxtbml(f'{self.wd}/{name}/ml', **self.xtb_setting)
        # update ml results to rec
        for name in reactants:
            if "ml_feature" not in rec[name].keys():
                try:
                    ml_feature = pd.read_csv(f'{self.wd}/{name}/ml/ml_feature.csv')
                    if 'HOAO_a (eV)' in ml_feature.columns:
                        ml_df = ml_feature[['atomic partial charges', 'coordination number',
                                            'HOAO_a (eV)', 'HOAO_b (eV)',
                                            'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]

                    else:
                        ml_df = ml_feature[['atomic partial charges', 'coordination number',
                                            'HOAO (eV)', 'LUAO (eV)', 'Atom']]
                    rec[name]['ml_feature'] = ml_df
                except:
                    sys.exit(f'{self.wd}/{name}/ml wrong')

                try:
                    if os.path.getsize(f'{self.wd}/{name}/ml/wboe') != 0:
                        wboe = pd.read_csv(f'{self.wd}/{name}/ml/wboe', header=None, sep='\s+')
                        rec[name]['wboe'] = wboe
                    else:
                        rec[name]['wboe'] = []
                except:
                    sys.exit(f'{self.wd}/{name}/wboe wrong')

        reactants_total_dict = {}
        # association
        # get target asoociation pairs
        association_d = {}
        association_d, hoaoa_df, hoaob_df, luaoa_df, luaob_df = association_pairs(reactants, rec,
                                                                                  **self.searching_setting)
        ninter = len(association_d.keys())
        # doscan
        reactor = self.searching_setting.get('reactor_like', False)
        reactants_inter_dict = push(association_d, rec, self.wd, reactor, stepnr, **self.xtb_setting)

        # dissociation
        # get target dissociation pairs
        dissociation_d = {}
        dissociation_d = dissociation_pair(reactants, rec, **self.searching_setting)

        # do scan
        reactants_dis_dict = pull(dissociation_d, rec, self.wd, reactor, ninter, stepnr, **self.xtb_setting)

        # # check scan results
        del_l = []
        reactants_total_dict = {**reactants_inter_dict, **reactants_dis_dict}
        if reactants_total_dict:
            path_l = list(reactants_total_dict.keys())
            for i in path_l:
                del_l += check_xtb_scan(self.wd, i)
            tmp = killpath(rec, del_l, reactants_total_dict)
            reactants_total_dict, rec = tmp[0],tmp[1]

        # subprocess: fragementation, assigning charge, m and get filed path due to gfn2 opt
        del_l = []
        if reactants_total_dict:
            path_l = list(reactants_total_dict.keys())
            for i in path_l:
                del_l_tmp, rec_tmp = subprocess(self.wd, i, reactants, reactants_total_dict, rec, reactor,
                                                **self.kwargs)
                del_l += del_l_tmp
                rec = update_dict(rec, rec_tmp)
            tmp = killpath(rec, del_l, reactants_total_dict)
            reactants_total_dict, rec = tmp[0], tmp[1]

        # get sp energy difference
        if reaction_rec is None:
            reaction_rec = {}
        if self.orca_setting.get("dft_level", 'gfn2') != 'gfn2':
            typ = 'dft'
        else:
            typ = 'gfn2'

        if reactants_total_dict:
            path_l = list(reactants_total_dict.keys())
            del_l = []
            call_header("SP", " ")
            for path_i in path_l:
                reaction_rec_tmp, del_l_tmp = elvauate_reaction(path_i, self.wd, rec, reactants_total_dict,
                                                                reaction_rec, typ, **self.searching_setting)

                reaction_rec = update_dict(reaction_rec, reaction_rec_tmp)
                del_l = del_l + del_l_tmp
            print(f"reaction_rec for step {stepnr}: {reaction_rec.keys()}")
            print(f"del_l for step {stepnr}: {del_l}")
        reactants_total_dict, rec,reaction_rec = killpath(rec, del_l, reactants_total_dict,reaction_rec)

        # dft opt
        if self.orca_setting.get("dft_opt", False) and reactants_total_dict:
            path_l = list(reactants_total_dict.keys())
            call_header("SPOPT", " ")
            for path_i in path_l:
                fragment_l = [f'{path_i}/{k}' for k in os.listdir(f'{self.wd}/{path_i}') if
                              os.path.isdir(f'{self.wd}/{path_i}/{k}') and str.isdigit(k)
                              and f'{path_i}/{k}' in rec.keys()]
                for struc in fragment_l:
                    m = int(rec[f'{struc}']['m'])
                    charge = rec[f'{struc}']['charge']
                    nat = rec[f'{struc}']['coord'].shape[0]
                    if nat > 1:
                        if m != 2:
                            rc = Struc(rec[f'{struc}']['coord'], m, charge). \
                                runorcasp(f'{self.wd}/{struc}/orca/orca_{str(m)}', caltyp='sp OPT', **self.orca_setting)
                            if rc != 0:
                                sys.exit(f'{self.wd}/{struc}/orca/orca_{str(m)} sp opt wrong')
                            rec = collect_sp_xyz(rec, f'{struc}', f'{self.wd}/{struc}/orca/orca_{str(m)}', 'dft')
                        else:
                            rc = Struc(rec[f'{struc}']['coord'], m, charge). \
                                runorcasp(f'{self.wd}/{struc}/orca', caltyp='sp OPT', **self.orca_setting)
                            if rc != 0:
                                sys.exit(f'{self.wd}/{struc}/orca sp opt wrong')
                            rec = collect_sp_xyz(rec, f'{struc}', f'{self.wd}/{struc}/orca/', 'dft')
                # collect reaction energy again
                reaction_rec_tmp, del_l_tmp = elvauate_reaction(path_i, self.wd, rec, reactants_total_dict,
                                                                reaction_rec, typ, **self.searching_setting)
                reaction_rec = update_dict(reaction_rec, reaction_rec_tmp)
        # do trv and solv
        if self.orca_setting.get("dft_opt", False) and (self.orca_setting.get("reaction_energy", 'sp') == 'g' or
            self.orca_setting.get("reaction_energy", 'sp') == 'G') and reactants_total_dict:
            path_l = list(reactants_total_dict.keys())
            call_header("G", " ")
            for path_i in path_l:
                fragment_l = [f'{path_i}/{k}' for k in os.listdir(f'{self.wd}/{path_i}') if
                              os.path.isdir(f'{self.wd}/{path_i}/{k}') and str.isdigit(k)
                              and f'{path_i}/{k}' in rec.keys()]
                for struc in fragment_l:
                    m = int(rec[f'{struc}']['m'])
                    charge = rec[f'{struc}']['charge']
                    rc1 = Struc(rec[f'{struc}']['coord'], m, charge). \
                        runtrv(f'{self.wd}/{struc}/gfn2_trv', **self.xtb_setting)
                    rc2 = Struc(rec[f'{struc}']['coord'], m, charge). \
                        runsolv(f'{self.wd}/{struc}/gfn2_solv', **self.xtb_setting)
                    rec = collect_g_xyz(rec, f'{struc}', f'{self.wd}/{struc}')
                # collect reaction energy again
                reaction_rec_tmp, del_l_tmp = elvauate_reaction(path_i, self.wd, rec, reactants_total_dict,
                                                                reaction_rec, 'g', **self.searching_setting)
                reaction_rec = update_dict(reaction_rec, reaction_rec_tmp)

        return reactants_total_dict, rec, reaction_rec, \
               (association_d, dissociation_d, hoaoa_df, hoaob_df, luaoa_df, luaob_df)


def killpath(rec, del_l, reactants_total_dict,reaction_rec=None):
    for i in list(set(del_l)):
        del reactants_total_dict[i]
        for k in list(rec.keys()):
            if k.split('/')[0] == i:
                del rec[k]
        if reaction_rec is not None:
            if i in reaction_rec.keys():
                del reaction_rec[i]
    # if reactor:
    #     reactant_l = reactants_total_dict[i]  # {'11': ['1/3', '1/5'], '12': ['1/3', '1/5'], '13': ['1/1']}
    #
    # else:
    #     reactant_l = [k for k in rec.keys() if
    #                   k.split('/')[0] == reactants_total_dict[i][0].split('/')[0]]
    return reactants_total_dict, rec,reaction_rec
