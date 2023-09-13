import builtins
import os
import shutil
import sys
import time
import copy
import traceback
from functools import partial

import billiard as multiprocessing
from Assign_reaction import do_assign_reactor, do_reaction
from operation import write_rec
from read import update_dict
from Subprocess import elvauate_reaction
from collections import Counter


def reactionnetwork(rec, wd, **kwargs):
    # rec: Record all reactants info: name, coord, at_features, sp
    # rec_valid: Record non-repeat reactants info: name, coord, at_features, sp
    # reaction_rec: Record valid reactions paths {path_name: [[R],[P],[delt_E]}
    # failreaction_rec: Record failed reactions paths {nr: [[R],[P]}
    # wd : Working folder
    # cycle: max simulation cycles

    if kwargs.get('kwargs'):
        kwargs = kwargs.get('kwargs')
    cycle = kwargs.get('cycle', 1)
    reactor = kwargs.get('reactor_like', False)
    energy_lim = kwargs.get('elim', 200)
    num_cores = kwargs.get('num_cores', 4)
    pid = kwargs.get('pid', None)
    if kwargs.get('dft_level'):
        typ = 'dft'
    else:
        typ = 'gfn2'
    print(f"main pid={pid}")

    reaction_rec = {}
    failreaction_rec = {}
    rec_valid = copy.deepcopy(rec)
    # if "refproduct" in rec_valid.keys():
    #     del rec_valid["refproduct"]
    step = 0
    while step < cycle:
        if reactor:
            react_l = [k for k in list(rec_valid.keys())]  # all species appeared in system
            print(f"nr of current species: {len(react_l)}")
            print(f"rec_l for step {step}: {react_l} ")
            t1 = time.time()
            reactants_tmp, fail_tmp = do_assign_reactor(react_l, rec, reaction_rec,
                                                        failreaction_rec)  # random select two molecules
            failreaction_rec = update_dict(failreaction_rec, fail_tmp)
            print(f"reactants = {reactants_tmp}")
            t2 = time.time()
            print(f'for cycle {step}, do_assign_reactor costs {(t2 - t1) / 60} min')
            reactants_total_dict, rec_tmp, fail_tmp = do_reaction(wd, rec, reactants_tmp, failreaction_rec, step + 1,
                                                                  kwargs=kwargs)  # 1 asso + 1 disso
            rec = update_dict(rec, rec_tmp)
            failreaction_rec = update_dict(failreaction_rec, fail_tmp)

        else:
            print(f"rec_key: {rec.keys}")
            reactants_total_dict = {}
            path_l = list(set([k.split('/')[0] for k in list(rec.keys())
                               if len(k.split('/')[0]) == step + 1
                               and os.path.isdir(f'{wd}/{k.split("/")[0]}')]))  # reaction path for this cycle
            print(f"path for this cycle:{path_l}")
            # individually sim each path
            items = []
            for i in range(len(path_l)):
                path_name = path_l[i]
                reactants_tmp = [k for k in list(rec.keys())if k.split("/")[0] == path_name]
                print(f"reactants for this path {path_name}: {reactants_tmp}")
                items.append((wd, rec, reactants_tmp, None,0))

            with multiprocessing.Pool(num_cores) as pool:
                mapfunc = partial(do_reaction, **kwargs)
                for result in pool.starmap(mapfunc, items):
                    reactants_total_dict = update_dict(reactants_total_dict, list(result)[0])
                    rec = update_dict(rec, list(result)[1])
                    failreaction_rec = update_dict(failreaction_rec, list(result)[2])
                    # print(f'Got result for reactants_total_dict: {reactants_total_dict}', flush=True)

        if reactants_total_dict:
            # collect energy and evaluate reaction
            path_l = list(reactants_total_dict.keys())
            del_l = []
            for path_i in path_l:
                reaction_rec_tmp, del_l_tmp = elvauate_reaction(path_i, wd, rec, energy_lim, reactants_total_dict,
                                                                reaction_rec, typ,reactor)
                reaction_rec = update_dict(reaction_rec, reaction_rec_tmp)
                del_l = del_l + del_l_tmp
            print(f"reaction_rec for step {step}: {reaction_rec.keys()}")
            print(f"del_l for step {step}: {del_l}")

            # kill paths that disqualified
            for i in list(set(del_l)):
                if os.path.isdir(f'{wd}/{i}') and reactor:
                    shutil.rmtree(f'{wd}/{i}', ignore_errors=True)
                for j in list(rec.keys()):
                    if j.split('/')[0] == i:
                        del rec[j]
                failreaction_rec[len(failreaction_rec.keys()) + 1] = reaction_rec[i]
                del reaction_rec[i]
                path_l.remove(i)

            # update rec_valid
            rec_valid = {}
            tmp = []
            for i in list(rec.keys()):
                if 'mb' in rec[i].keys():
                    mb = rec[i]['mb']
                    if mb not in tmp:
                        tmp.append(mb)
                        rec_valid[i] = {}
                        rec_valid[i] = rec[i]
            # write document
            if reactor:
                if path_l:  # step +1 only if there is association [s for s in path_l if "_1" in s]
                    write_rec(rec_valid, wd, 'molecules')
                    write_rec(reaction_rec, wd, 'reaction')
                    write_rec(failreaction_rec, wd, 'failreaction')
                    new_n_reactants = len(list(rec_valid.keys()))
                    sussnr = len([k for k in reaction_rec.keys() if '_1' in k])
                    failnr = len([k for k in failreaction_rec.keys() if len(failreaction_rec[k][0]) > 1])
                    if failnr >= new_n_reactants ** 2 - sussnr:
                        print("no more possible reactions")
                        break
                    else:
                        step += 1
                else:
                    continue
            else:
                if path_l:
                    # write_rec(rec_valid, wd, 'molecules')
                    write_rec(rec, wd, 'rec')
                    write_rec(reaction_rec, wd, 'reaction')
                    write_rec(failreaction_rec, wd, 'failreaction')
                    # tmp test
                    product_mb = [rec[k]['mb'] for k in rec.keys() if 'refproduct' in k]

                    for k in list(reaction_rec.keys()):
                        products_k = reaction_rec[k][1]
                        if Counter(product_mb) == Counter(products_k):
                            print("!!!!!!!found product!!!!!!!!")
                            try:
                                if os.path.isdir(f'{wd}/{k}'):
                                    shutil.move(f'{wd}/{k}',f'{wd}/{k}_hier')
                            except:
                                traceback.print_exc()
                                pass
                            break
                    print("iamhier", step)
                    step += 1

                else:
                    print("no path_l")
                    break
    else:
        # write_rec(rec, wd, 'rec')
        # write_rec(reaction_rec, wd, 'reaction')
        # write_rec(failreaction_rec, wd, 'failreaction')
        # product_mb = [rec[k]['mb'] for k in rec.keys() if 'refproduct' in k]
        # products_d = {}
        # for k in list(reaction_rec.keys()):
        #     products_d[reaction_rec[k][1]] = reaction_rec[k][3]
        # products_l = [reaction_rec[k][1] for k in list(reaction_rec.keys())]
        # for k in products_l:
        #     if Counter(product_mb) == Counter(k):
        #         # if product_mb in products_l:
        #         print("!!!!!!!found product!!!!!!!!")
        #         try:
        #             path_i = products_d[k]
        #             if os.path.isdir(f'{wd}/{path_i}'):
        #                 shutil.move(f'{wd}/{path_i}', f'{wd}/{path_i}_hier')
        #         except:
        #             traceback.print_exc()
        #             pass

        print("done")
        return rec
