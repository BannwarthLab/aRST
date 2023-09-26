import builtins
import os
import shutil

import time
import copy
import traceback
from functools import partial
import billiard as multiprocessing
from multiprocessing import set_start_method
from Assign_reaction import do_assign_reactor,Singlestep
from operation import write_rec
from read import update_dict
from collections import Counter


# def reactionnetwork(rec, wd, reaction_rec=None, failreaction_rec=None, step=0, **kwargs):
#     # rec: Record all reactants info: name, coord, at_features, sp
#     # rec_valid: Record non-repeat reactants info: name, coord, at_features, sp
#     # reaction_rec: Record valid reactions paths {path_name: [[R],[P],[delt_E]}
#     # failreaction_rec: Record failed reactions paths {nr: [[R],[P]}
#     # wd : Working folder
#     # cycle: max simulation cycles
#     searching_setting = kwargs.get('searching_setting')
#     cycle = searching_setting.get('cycle', 1)
#     reactor = searching_setting.get('reactor_like', False)
#     num_cores = searching_setting.get('num_cores', 4)
#
#     def do_reaction_single(reactants, step, rec, wd, kwargs):
#         singlestep_instance = Singlestep(wd, **kwargs)
#         return singlestep_instance.do_reaction(reactants, step, rec)
#     reactants = ["1/1", "1/2"]
#     wd = "1"
#     step = 1
#     kwargs = {"searching_setting": {"orca_setting": 1, "xtb_setting": 2}}
#     num_cores = 1
#     rec = {}
#     items = [["1/1", "1/2"], ["1/1", "1/2"], ["1/1"]]
#     import billiard as multiprocessing
#     with multiprocessing.Pool(num_cores) as pool:
#         do_reaction_partial = partial(do_reaction_single, wd=wd, kwargs=kwargs)
#         results = pool.starmap(do_reaction_partial,
#                                [(reactants_tmp, step, rec) for reactants_tmp in items])
#     print(results)
#
#     # while step < cycle:
#     #     if reactor:
#     #         react_l = [k for k in list(rec.keys())]  # all species appeared in system
#     #         print(f"nr of current species: {len(react_l)}")
#     #         print(f"rec_l for step {step}: {react_l} ")
#     #         t1 = time.time()
#     #         reactants_tmp, fail_tmp = do_assign_reactor(react_l, rec, reaction_rec,
#     #                                                     failreaction_rec)  # random select two molecules
#     #         failreaction_rec = update_dict(failreaction_rec, fail_tmp)
#     #         print(f"reactants = {reactants_tmp}")
#     #         t2 = time.time()
#     #         print(f'for cycle {step}, do_assign_reactor costs {(t2 - t1) / 60} min')
#     #
#     #         reactants_total_dict_tmp, rec_tmp, reaction_rec_tmp = \
#     #             Singlestep(wd, **kwargs).do_reaction(reactants_tmp,step, rec, reaction_rec)
#     #
#     #         rec = update_dict(rec, rec_tmp)
#     #         reaction_rec = update_dict(reaction_rec, reaction_rec_tmp)
#     #
#     #
#     #     else:
#     #         print(f"rec_key: {rec.keys}")
#     #         reactants_total_dict = {}
#     #         path_l = list(set([k.split('/')[0] for k in list(rec.keys())
#     #                            if len(k.split('/')[0]) == step + 1
#     #                            and os.path.isdir(f'{wd}/{k.split("/")[0]}')]))  # reaction path for this cycle
#     #         print(f"path for this cycle:{path_l}")
#     #
#     #         # for i in range(len(path_l)):
#     #         #     path_name = path_l[i]
#     #         #     reactants_tmp = [k for k in list(rec.keys()) if k.split("/")[0] == path_name]
#     #         #     reactants_total_dict_tmp, rec_tmp, reaction_rec_tmp = \
#     #         #         Singlestep(wd, **kwargs).do_reaction(reactants_tmp, step, rec, reaction_rec)
#     #         #     reactants_total_dict = update_dict(reactants_total_dict, reactants_total_dict_tmp)
#     #         #     rec = update_dict(rec, rec_tmp)
#     #         #     reaction_rec = update_dict(reaction_rec, reaction_rec_tmp)
#     #         # path_l = list(set([k.split('/')[0] for k in list(rec.keys())
#     #         #                    if len(k.split('/')[0]) == step + 1
#     #         #                    and os.path.isdir(f'{wd}/{k.split("/")[0]}')]))  # reaction path surviored for this cycle
#     #
#     #         # individually sim each path
#     #         items = []
#     #         for i in range(len(path_l)):
#     #             path_name = path_l[i]
#     #             reactants_tmp = [k for k in list(rec.keys()) if k.split("/")[0] == path_name]
#     #             print(f"reactants for this path {path_name}: {reactants_tmp}")
#     #             items.append(reactants_tmp)
#     #
#     #         def do_reaction_single(reactants, step, rec, wd, kwargs):
#     #             singlestep_instance = Singlestep(wd, **kwargs)
#     #             return singlestep_instance.do_reaction(reactants, step, rec)
#     #         print("iamhier!!!!@@@#")
#     #         reactants = ["1/1", "1/2"]
#     #         wd = "1"
#     #         step = 1
#     #         kwargs = {"searching_setting": {"orca_setting": 1, "xtb_setting": 2}}
#     #         num_cores = 1
#     #         rec = {}
#     #         items = [["1/1", "1/2"], ["1/1", "1/2"], ["1/1"]]
#     #         # Use multiprocessing.Pool to parallelize the tasks
#     #         with multiprocessing.Pool(1) as pool:
#     #             do_reaction_partial = partial(do_reaction_single, wd=wd, kwargs=kwargs)
#     #             results = pool.starmap(do_reaction_partial,
#     #                                    [(reactants_tmp, step, rec) for reactants_tmp in items])
#     #         print(len(results))
#     #         #     for result in results:
#     #         #         reactants_total_dict = update_dict(reactants_total_dict, list(result)[0])
#     #         #         rec = update_dict(rec, list(result)[1])
#     #         #         reaction_rec = update_dict(reaction_rec, list(result)[2])
#     #         # pool.close()
#     #         # pool.join()
#     #         path_l = list(set([k.split('/')[0] for k in list(rec.keys())
#     #                            if len(k.split('/')[0]) == step + 1
#     #                            and os.path.isdir(f'{wd}/{k.split("/")[0]}')]))  # reaction path surviored for this cycle
#     #         # print(f'Got result for reactants_total_dict: {reactants_total_dict}', flush=True)
#     #
#     #         # write document
#     #         # if reactor:
#     #         #     if path_l:  # step +1 only if there is association [s for s in path_l if "_1" in s]
#     #         #         write_rec(reaction_rec, wd, 'reaction')
#     #         #         write_rec(failreaction_rec, wd, 'failreaction')
#     #         #         sussnr = len([k for k in reaction_rec.keys() if '_1' in k])
#     #         #         # failnr = len([k for k in failreaction_rec.keys() if len(failreaction_rec[k][0]) > 1])
#     #         #         # if failnr >= new_n_reactants ** 2 - sussnr:
#     #         #         #     print("no more possible reactions")
#     #         #         #     break
#     #         #         # else:
#     #         #         if sussnr:
#     #         #             step += 1
#     #         #     else:
#     #         #         continue
#     #         # else:
#     #         if path_l:
#     #             write_rec(rec, wd, 'rec')
#     #             write_rec(reaction_rec, wd, 'reaction')
#     #             # write_rec(failreaction_rec, wd, 'failreaction')
#     #             # tmp test
#     #             product_mb = [rec[k]['mb'] for k in rec.keys() if 'refproduct' in k]
#     #
#     #             for k in list(reaction_rec.keys()):
#     #                 products_k = reaction_rec[k][1]
#     #                 if Counter(product_mb) == Counter(products_k):
#     #                     print("!!!!!!!found product!!!!!!!!")
#     #                     try:
#     #                         if os.path.isdir(f'{wd}/{k}'):
#     #                             shutil.move(f'{wd}/{k}', f'{wd}/{k}_hier')
#     #                     except:
#     #                         traceback.print_exc()
#     #                         pass
#     #                     break
#     #             print("iamhier", step)
#     #             step += 1
#     #
#     #         else:
#     #             print("no path_l")
#     #             break
#     # else:
#     #     write_rec(rec, wd, 'rec')
#     #     write_rec(reaction_rec, wd, 'reaction')
#     #     write_rec(failreaction_rec, wd, 'failreaction')
#     #     product_mb = [rec[k]['mb'] for k in rec.keys() if 'refproduct' in k]
#     #     products_d = {}
#     #     for k in list(reaction_rec.keys()):
#     #         products_d[reaction_rec[k][1]] = reaction_rec[k][3]
#     #     products_l = [reaction_rec[k][1] for k in list(reaction_rec.keys())]
#     #     for k in products_l:
#     #         if Counter(product_mb) == Counter(k):
#     #             # if product_mb in products_l:
#     #             print("!!!!!!!found product!!!!!!!!")
#     #             try:
#     #                 path_i = products_d[k]
#     #                 if os.path.isdir(f'{wd}/{path_i}'):
#     #                     shutil.move(f'{wd}/{path_i}', f'{wd}/{path_i}_hier')
#     #             except:
#     #                 traceback.print_exc()
#     #                 pass
#     #
#     #     print("done")
#     #     return rec


