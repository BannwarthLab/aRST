import time
from Subprocess import subprocess, collect_dft
from read import existreaction, update_dict
import random
from PushAndPull import push, pull
from run_mlcalculation import ml, collect_ml


def do_assign_reactor(react_l, rec, reaction_rec, fail_rec):
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
            # test if reaction already exist
            if existreaction(fail_rec, reactants_mb):
                reactants_tmp = []
            else:
                if existreaction(reaction_rec, reactants_mb):
                    reactants_tmp = []
                else:
                    # test if unusual reactants
                    if tot_c > 1 or tot_c < -1:
                        fail_rec[len(fail_rec.keys()) + 1] = reactants_mb
                        reactants_tmp = []
    else:
        reactants_tmp = []
    return reactants_tmp, fail_rec


def do_reaction(wd, rec, react_l, fail_rec=None, step=0, **kwargs):
    if kwargs.get('kwargs'):
        kwargs = kwargs.get('kwargs')
    xtbwd = kwargs.get('xtbwd')
    alpb = kwargs.get('alpb', None)
    reactor = kwargs.get('reactor_like', False)
    level = kwargs.get('dft_level', None)
    mode = kwargs.get('platform', 0)

    t1 = time.time()
    # ml features
    ml(rec, react_l, wd, xtbwd, alpb)
    rec = collect_ml(rec, react_l, wd)

    # Asoociation
    if len(react_l) > 1:
        reactants_inter_dict = push(react_l, rec, wd, step=step, kwargs=kwargs)
        ninter = len(reactants_inter_dict.keys())
    else:
        reactants_inter_dict = {}
        ninter = 0

    # Dissociation
    kwargs.update({'association_count': ninter})
    reactants_dis_dict = pull(react_l, rec, wd, step=step, kwargs=kwargs)
    del kwargs['association_count']
    t2 = time.time()

    # check scan stauts
    reactants_inter_dict.update(reactants_dis_dict)  # now inter_dict = total dict
    if reactants_inter_dict:
        path_l = list(reactants_inter_dict.keys())
        del_l = []
        for i in path_l:
            rec_tmp,del_l_tmp = subprocess(wd, i, react_l, reactants_inter_dict, rec, level, alpb, reactor, mode)

            rec = update_dict(rec, rec_tmp)
            del_l += del_l_tmp

        for i in list(set(del_l)):
            path_l.remove(i)
            if reactor:
                reactant_l = reactants_inter_dict[i]  # {'11': ['1/3', '1/5'], '12': ['1/3', '1/5'], '13': ['1/1']}

            else:
                reactant_l = [k for k in rec.keys() if k.split('/')[0] == reactants_inter_dict[i][0].split('/')[0]]

            # if fail_rec is not None:
            #     fail_rec[len(fail_rec.keys()) + 1] = [[rec[k]['mb'] for k in reactant_l ],0,0]
            # else:
            #     fail_rec = {}
            #     fail_rec[1] = [[rec[k]['mb'] for k in reactant_l], 0, 0]
            del reactants_inter_dict[i]
            for k in list(rec.keys()):
                if k.split('/')[0] == i:
                    del rec[k]

        if path_l and level is not None:
            rec = update_dict(rec, collect_dft(rec, wd))

    t3 = time.time()
    return reactants_inter_dict, rec, fail_rec  # inter_dict = total dict
