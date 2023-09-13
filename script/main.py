import builtins
import os
import shutil
from argparse import ArgumentParser
from read import read_input, prepare, update_dict
from cycle import reactionnetwork
import time
import billiard as multiprocessing
from multiprocessing import set_start_method
from functools import partial
import sys


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


    kwargs, rec = read_input(args.in_file)
    print(kwargs)
    if kwargs.get('platform') == 0:
        num_cores = int(os.getenv('OMP_NUM_THREADS'))
        kwargs['num_cores'] = num_cores

    else:
        num_cores = 4
        kwargs['num_cores'] = 4
    if not rec:
        print("Failed to read input")
        sys.exit()

    # determine the working folder
    wd = os.path.abspath(os.path.dirname(args.in_file))
    print(f"working folder is {wd}")
    print("num_cores: ",num_cores)

    print(sys.version)
    print(sys.path)

    # clean up old resultes
    l = os.listdir(wd)
    for i in l:
        if os.path.isdir(f'{wd}/{i}'):
            shutil.rmtree(f'{wd}/{i}')
        elif i.endswith('out'):
            os.remove(f'{wd}/{i}')
        else:
            if i.endswith('json'):
                os.remove(f'{wd}/{i}')

    # recording starting time
    t1 = time.time()

    # prepare reactants before cycle: do opt and sp and get basic propoties for each reactant
    with multiprocessing.Pool(num_cores) as pool:
        items = [(wd, i, rec) for i in list(rec.keys())]
        mapfunc = partial(prepare, **kwargs)
        for result in pool.starmap(mapfunc, items):
            rec = update_dict(rec, result)

    # begin cycle
    pid = str(os.getpid())
    print(pid)
    kwargs['pid'] = pid
    rec = reactionnetwork(rec, wd, kwargs=kwargs)

    # recording finishing time
    t2 = time.time()
    print(f"finished with {t2 - t1}")
    sys.exit()
