import os
import shutil
import sys
import time
import traceback
import billiard as multiprocessing
from argparse import ArgumentParser
from read import read_input, prepare, update_dict
from multiprocessing import set_start_method
from functools import partial
from head import print_header, call_header, printout_cycle_results
from Assign_reaction import do_assign_reactor, Singlestep
from operation import write_rec
from collections import Counter

if __name__ == "__main__":
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

    import warnings

    warnings.filterwarnings("ignore")
    set_start_method('spawn')

    # get input
    rec, kwargs = read_input(args.in_file)

    # set cores
    num_cores = int(os.getenv('OMP_NUM_THREADS', '1'))

    if kwargs.get("searching_setting").get("platform") == 0:
        num_cores = 4
    kwargs['searching_setting']['num_cores'] = num_cores

    if not rec:
        print("Failed to read input")
        sys.exit()

    # determine the working folder
    wd = os.path.abspath(os.path.dirname(args.in_file))
    print(f"working folder is {wd}")
    print("num_cores: ", kwargs['searching_setting']['num_cores'])

    # print(sys.version)
    # print(sys.path)

    # prepare reactants before cycle: do opt and sp and get basic propoties for each reactant
    items = [(wd, i, rec) for i in list(rec.keys()) if "sp_gfn2" not in rec[i].keys()]
    if 0 < len(items) < num_cores:
        nmp = len(items)
    else:
        nmp = num_cores
    with multiprocessing.Pool(nmp) as pool:
        mapfunc = partial(prepare, **kwargs)
        for result in pool.starmap(mapfunc, items):
            rec = update_dict(rec, result)

    # recording starting time
    t1 = time.time()
    pid = str(os.getpid())
    # print("pid:", pid)

    # begin cycle
    del kwargs['geom']
    print_header()

    # rec = reactionnetwork(rec, wd, **kwargs)
    searching_setting = kwargs.get('searching_setting')
    cycle = searching_setting.get('cycle', 1)
    step = 0
    reaction_rec = {}
    print('Reactants: ' + ''.join([rec[i]['mb'].split('|')[2] + str(" ; ")
                                   for i in list(rec.keys()) if i.split('/')[0] == '1']))
    while step < cycle:
        printout_d = {}
        call_header(f"cycle {step}", "=")
        reactants_total_dict = {}
        path_l = list(set([k.split('/')[0] for k in list(rec.keys())
                           if len(k.split('/')[0]) == step + 1
                           and os.path.isdir(f'{wd}/{k.split("/")[0]}')]))  # reaction path for this cycle
        print(f"path for this cycle:{path_l}")

        # individually sim each path
        items = []
        for i in range(len(path_l)):
            path_name = path_l[i]
            reactants_tmp = [k for k in list(rec.keys()) if k.split("/")[0] == path_name]
            items.append(reactants_tmp)
            printout_d[path_name] = {}
            printout_d[path_name]['reactants'] = reactants_tmp


        def do_reaction_single(reactants, step, rec, wd, kwargs):
            singlestep_instance = Singlestep(wd, **kwargs)
            return singlestep_instance.do_reaction(reactants, step, rec)


        # Use multiprocessing.Pool to parallelize the tasks
        if 0 < len(items) < num_cores:
            nmp = len(items)
        else:
            nmp = num_cores
        with multiprocessing.Pool(nmp) as pool:
            do_reaction_partial = partial(do_reaction_single, wd=wd, kwargs=kwargs)
            results = pool.starmap(do_reaction_partial,
                                   [(reactants_tmp, step, rec) for reactants_tmp in items])
        for r in range(len(list(results))):
            result = list(results)[r]
            reactants_total_dict = update_dict(reactants_total_dict, result[0])
            rec = update_dict(rec, result[1])
            reaction_rec = update_dict(reaction_rec, result[2])
            # for print out
            if reactants_total_dict:
                reactants_tmp = reactants_total_dict[list(reactants_total_dict.keys())[0]]
                path_i = list(reactants_total_dict.keys())[0].strip()
                association_d, hoaoa_df, hoaob_df, luaoa_df, luaob_df = result[3]
                printout_cycle_results(path_i, reactants_tmp, association_d, hoaoa_df, hoaob_df, luaoa_df,
                                       luaob_df, rec)

        path_l = list(set([k.split('/')[0] for k in list(rec.keys())
                           if len(k.split('/')[0]) == step + 2
                           and os.path.isdir(f'{wd}/{k.split("/")[0]}')]))  # reaction path surviored for this cycle

        if path_l:
            write_rec(rec, wd, 'rec')
            write_rec(reaction_rec, wd, 'reaction')
            # write_rec(failreaction_rec, wd, 'failreaction')
            # tmp test
            product_mb = [rec[k]['mb'] for k in rec.keys() if 'refproduct' in k]

            for k in list(reaction_rec.keys()):
                products_k = reaction_rec[k][1]
                if Counter(product_mb) == Counter(products_k):
                    print("!!!!!!!found product!!!!!!!!")
                    try:
                        if os.path.isdir(f'{wd}/{k}'):
                            shutil.move(f'{wd}/{k}', f'{wd}/{k}_hier')
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
        print("step=", step)
        write_rec(rec, wd, 'rec')
        write_rec(reaction_rec, wd, 'reaction')

        print("done")

    # recording finishing time
    t2 = time.time()
    print(f"finished with {t2 - t1}")
    sys.exit()
