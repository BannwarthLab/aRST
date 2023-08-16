from Subprocess import subprocess, collect_dft
from read import existreaction, update_dict
import random
import billiard as multiprocessing
from PushAndPull import push, pull
from run_mlcalculation import ml, collect_ml


def do_assign_reactor(react_l, rec, reaction_rec):
    reactants_tmp = []  # reactants for this cycle
    if len(react_l) > 1:
        while not reactants_tmp:  # assume that +2/-2 charged molecules are not stable and can only dissociate back
            reactants_tmp = random.sample(react_l, 2)
            reactants_mb = [rec[k]['mb'] for k in reactants_tmp]
            tot_c = float('inf')
            try:
                tot_c = sum([rec[k]['charge'] for k in reactants_tmp])
            except:
                pass
            # test if reaction already exist
            if existreaction(reaction_rec, reactants_mb):
                reactants_tmp = []
            else:
                # test if unusual reactants
                if tot_c > 1 or tot_c < -1:
                    reactants_tmp = []
    else:
        reactants_tmp = react_l
    return reactants_tmp


def do_reaction(wd, rec, react_l, step=0, **kwargs):
    if kwargs.get('kwargs'):
        kwargs = kwargs.get('kwargs')

    xtbwd = kwargs.get('xtbwd')
    alpb = kwargs.get('alpb', None)
    reactor = kwargs.get('reactor_like', False)
    level = kwargs.get('dft_level', 'r2SCAN-3c')

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

    # check scan stauts
    reactants_inter_dict.update(reactants_dis_dict)  # now inter_dict = total dict
    path_l = list(reactants_inter_dict.keys())
    with multiprocessing.Pool(len(path_l)) as pool:
        # pool = multiprocessing.Pool(len(path_l))
        results = []
        for i in path_l:
            results.append(pool.apply_async(subprocess,
                                            (wd, i, react_l, reactants_inter_dict, rec, level, alpb, reactor)))
        pool.close()
        pool.join()

    del_l = []  # list for path where scan failed
    for res in results:
        rec = update_dict(rec, list(res.get())[0])
        del_l += list(res.get())[1]

    for i in list(set(del_l)):
        path_l.remove(i)
        del reactants_inter_dict[i]

    pool = multiprocessing.Pool(len(path_l))
    results = [pool.apply_async(collect_dft, (rec, wd))]
    pool.close()
    pool.join()
    for res in results:
        rec = update_dict(rec, res.get())
    return reactants_inter_dict, rec  # inter_dict = total dict
