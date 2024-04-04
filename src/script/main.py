#!/usr/bin/python3.9
import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

import os
import shutil
import sys
import time
import traceback
# import billiard as multiprocessing
# from billiard import managers
from argparse import ArgumentParser
from read import read_input, prepare
from multiprocessing import set_start_method,log_to_stderr,get_logger,Manager
import multiprocessing

from functools import partial
from head import print_header, call_header, printout_cycle_results
from Assign_reaction import Singlestep,complete_struc_dftopt,complete_trvandsolv
from collections import Counter
import threading
from dictionaryInspector import DictionaryInspector,DataFrameEncoder

# Define a lock to make the function call thread-safe
function_lock = threading.Lock()


def do_reaction_single(reactants, reactionpath, shared_nstruc, shared_nreaction, STRUC, REACTION, wd, kwargs):
    with function_lock:
        singlestep_instance = Singlestep(wd, STRUC, REACTION, **kwargs)
        result = singlestep_instance.do_reaction(reactants, reactionpath, shared_nstruc, shared_nreaction)

        return result

def do_dftopt(wd,info_d,optlevel,pcm):
    with function_lock:
        result = complete_struc_dftopt(wd,info_d,optlevel,pcm)
        return result


if __name__ == "__main__":

    set_start_method('fork')

    # import logging

    # log_to_stderr()
    # get_logger().setLevel(logging.DEBUG)


    parser = ArgumentParser(
        # name of program
        prog="MEOW",
        # description
        description="mode 1: module load excstda/ orca/ xtb/; source ~/mb/env/bin/activate",
        # bottom text
        epilog="By Ying"
    )

    # add arguments

    parser.add_argument("in_file", type=str, help="The input file")

    args = parser.parse_args()



    # recording starting time
    t1 = time.time()

    # get input
    struc_rec, kwargs = read_input(args.in_file)
    for i in kwargs.keys():
        print(f"{i},{kwargs[i]}\n")
    if not struc_rec:
        print("Failed to read input")
        sys.exit()

    # set cores
    num_cores = int(os.getenv('OMP_NUM_THREADS', '1'))
    if kwargs.get("searching_setting").get("platform") == 0:
        num_cores = 4
    kwargs['searching_setting']['num_cores'] = num_cores

    # merge "searching_setting" and "xtb_stting"
    tmp_d = {**kwargs.get("searching_setting"), **kwargs.get("xtb_setting")}
    kwargs['searching_setting'] = tmp_d
    kwargs['xtb_setting'] = tmp_d

    # determine the working folder
    wd = os.path.abspath(os.path.dirname(args.in_file))
    print(f"working folder is {wd}")

    # print(sys.version)
    # print(sys.path)

    # prepare reactants before cycle: do opt and sp and get basic propoties for each reactant
    STRUC_inspector = DictionaryInspector(struc_rec)
    items = [(wd, struc_id, STRUC_inspector)
             for struc_id in struc_rec.keys() if
             STRUC_inspector.get_value(f"{struc_id}.struc_info.energy.sp_gfn2") is None]
    with multiprocessing.Pool(num_cores) as pool:
        mapfunc = partial(prepare, **kwargs)
        for result in pool.starmap(mapfunc, items):
            STRUC_inspector.update_strucrec_values(result)

    ini_reactants = [strucid for strucid in struc_rec.keys() if "ref" not in strucid]



    pid = str(os.getpid())
    # print("pid:", pid)

    # begin cycle
    del kwargs['geom']
    print_header()

    cycle = kwargs.get('searching_setting').get('cycle', 1)
    elim = kwargs.get('searching_setting').get('elim', 300)
    # elim evaluation
    if kwargs.get('orca_setting').get('reaction_energy', 'sp') == 'G':
        elim_mode = 'G'
    elif kwargs.get('orca_setting').get('dft_level', 'r2scan-3c') == 'gfn2':
        elim_mode = 'sp_gfn2'
    else:
        elim_mode = 'sp_dft'

    step = 0

    reaction_inspector = DictionaryInspector({})  # record all apeared reaction
    # set ref reaction, if there is one
    if "refp0001" in struc_rec.keys():
        initial_reactants = [i for i in struc_rec.keys() if "refp" not in i]
        initial_products = [i for i in struc_rec.keys() if "refp" in i]
        reactionid = "ref"
        reaction_inspector.set_reaction_values(STRUC_inspector, initial_reactants, initial_products,
                                               [], reactionid, reactionid,None,None)


    nstruc = len([k for k in struc_rec.keys() if not "ref" in k])
    nreaction = len([k for k in reaction_inspector.data.keys() if not "ref" in k])

    while step < cycle:

        # printout_d = {} not use for now
        call_header(f"cycle {step}", "=")
        # find paths that surviored last cycle and the corespongding reactants id
        items = [] # reacord reactants to be investigated for this step
        items2 = [] # reacord reaction path to be investigated for this step
        if step == 0:
            reactants = STRUC_inspector.find_sturcid_with_reactionpath('1')
            items.append(reactants)
            items2.append('1')
            # printout_d['1'] = {}
            # printout_d['1']['reactants'] = [STRUC_inspector.get_value(f'{i}.struc_info.general.formula')
            #                                 for i in reactants]
        else:
            rl = [] # reactionID list to be investigated for this step

            # start to screen reactionID
            for reactionid in reaction_inspector.data.keys():
                if reactionid != "ref":
                    reaction_pathinfo = reaction_inspector.get_value(f'{reactionid}.reaction_path')
                    # check if last step
                    lastcycle = False
                    path = None
                    if len(reaction_pathinfo) == 1:  # unique path
                        if len(reaction_pathinfo[0]) == step + 1:
                            lastcycle = True
                            path = reaction_pathinfo[0]
                    elif len(set([len(k) for k in reaction_pathinfo]))==1 and len(reaction_pathinfo[0]) == step + 1:
                        # degenerated reaction path in the same step
                        lastcycle = True
                        path = reaction_pathinfo[0]


                    if lastcycle:
                        # check reaction energy
                        ereaction = reaction_inspector.get_value(f'{reactionid}.energy.{elim_mode}')
                        if ereaction is None or ereaction<=elim:
                        # if ereaction < elim:
                            # check if new reactions
                            # i.e products_l + unreactants_l of this reactionID can't be any old reactants
                             if not reaction_inspector.test_if_reaction_loop(reactionid):
                                products_l = reaction_inspector.get_value(f"{reactionid}.struc_info.products") or []
                                unreactants_l = reaction_inspector.get_value(f"{reactionid}.struc_info.unreactants") or []

                                reactants = products_l + unreactants_l
                                items.append(reactants)
                                items2.append(path)

        if items:
            # Use Manager to creat two shared variante
            with Manager() as manager:
                #  int for counting struc/reaction number
                current_nstruc = manager.Value('i', nstruc)
                current_nreaction = manager.Value('i', nreaction)

                # Use multiprocessing.Pool to parallelize the tasks
                with multiprocessing.Pool(num_cores) as pool:
                    do_reaction_partial = partial(do_reaction_single, STRUC=STRUC_inspector,
                                                  REACTION=reaction_inspector, wd=wd, kwargs=kwargs)
                    results = pool.starmap(do_reaction_partial,
                                           [(reactants_tmp, path_tmp, current_nstruc, current_nreaction)
                                            for reactants_tmp, path_tmp in zip(items, items2)])
                    pool.close()
                    pool.join()


                    #### after reaction simulation

                # update shared values
                nstruc = current_nstruc.value
                nreaction = current_nreaction.value

            # receive results
            tmp_struc_d, tmp_reaction_d, out_d, out2_d = {}, {}, {}, {}
            for result in results:
                tmp_struc_d.update(result[0])
                tmp_reaction_d.update(result[1])
                out_d.update(result[2])
                out2_d.update(result[3])

            # kill reaction repeat within cycle (same reactants, different path, same products) and to old rections
            tmp_d = STRUC_inspector.update_struc_d_after_mp(tmp_struc_d)
            reaction_inspector.update_reaction_d_after_mp(tmp_reaction_d, tmp_d)


            try:
                # do ml for the new created struc
                Singlestep(wd,STRUC_inspector,reaction_inspector,**kwargs).ML_prepare(STRUC_inspector.data.keys())
            except:
                traceback.print_exc()


            toopt_struc_l = []

            for reactionid in out_d.keys():
                products = reaction_inspector.get_value(f"{reactionid}.struc_info.products") or []
                toopt_struc_l.extend(products)

            toopt_struc_l = [strucid for strucid in toopt_struc_l if strucid in STRUC_inspector.data.keys()]
            toopt_struc_l = list(set(toopt_struc_l))


            # do dft opt
            if kwargs.get('orca_setting').get('dft_opt', False):

                optlevel = kwargs.get("orca_setting", None).get("dft_level", None)
                pcm = kwargs.get("orca_setting", None).get("pcm", None)

                with multiprocessing.Pool(num_cores) as pool:
                    results = pool.starmap(do_dftopt,[(wd,{strucid: STRUC_inspector.data[strucid]},optlevel,pcm)
                                                      for strucid in toopt_struc_l])
                for result in results:
                    strucid = result[0]
                    STRUC_inspector.set_value(f"{strucid}.struc_info.general.coord",result[1])
                    STRUC_inspector.set_value(f"{strucid}.struc_info.general.multiplicity", result[2])
                    STRUC_inspector.set_value(f"{strucid}.struc_info.energy.sp_dft", result[3])

            # do trv and solv
            if kwargs.get('orca_setting').get('reaction_energy', 'sp') == 'G':
                xtb_setting = kwargs.get("xtb_setting")
                level = xtb_setting.get('xtblevel', 'gfn2')
                alpb = xtb_setting.get('alpb', 'water')
                with multiprocessing.Pool(num_cores) as pool:
                    results = pool.starmap(complete_trvandsolv,[(wd,{strucid: STRUC_inspector.data[strucid]},level,alpb)
                                                  for strucid in toopt_struc_l])
                for result in results:
                    strucid = result[0]
                    trv = result[1] if result[1] is not None else 0
                    solv = result[2] if result[2] is not None else 0
                    if kwargs.get('orca_setting').get('dft_level', 'gfn2') == 'gfn2':
                        sp = STRUC_inspector.get_value(f"{strucid}.struc_info.energy.sp_gfn2")
                        if sp is None:
                            sp =0
                    else:
                        sp = STRUC_inspector.get_value(f"{strucid}.struc_info.energy.sp_dft")
                        if sp is None:
                            sp = 0

                    STRUC_inspector.set_value(f"{strucid}.struc_info.energy.trv", trv)
                    STRUC_inspector.set_value(f"{strucid}.struc_info.energy.solv", solv)

                    STRUC_inspector.set_value(f"{strucid}.struc_info.energy.G",trv+solv+sp)


            # update reaction energy
            for reactionid,item in out_d.items():
                reactants_l = item[0]
                products_l = item[1]
                full_reactants = out2_d[reactionid]
                unreactants_l = [k for k in full_reactants if k not in reactants_l]
                if reactionid in reaction_inspector.data.keys():
                    reaction_inspector.set_reaction_values(STRUC_inspector,reactants_l,
                                                           products_l,None,reactionid,None,None,None,onlyenergy=True)


            step += 1
        else:

            step += 1
            # break

    # recording finishing time
    from toolbox import rm_folderORfile

    for folder in os.listdir(wd):
        if os.path.isdir(os.path.join(wd,folder)) and not ("struc" in folder or "scan" in folder):
            rm_folderORfile(os.path.join(wd,folder))

    import json

    # orgnize outputdata
    from dictionaryInspector import convert_df2json

    with open(f'{wd}/reaction.json', "w") as f:
        f.write(json.dumps(convert_df2json(reaction_inspector.data), indent=2))
        f.write('\n')

    with open(f'{wd}/struc.json', "w") as f:
        f.write(json.dumps(convert_df2json(STRUC_inspector.data), indent=2))
        f.write('\n')


    t2 = time.time()
    print(f"finished with {t2 - t1}")
    sys.exit()