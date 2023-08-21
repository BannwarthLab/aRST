import os
import shutil
import sys
from argparse import ArgumentParser
from read import read_input, prepare, update_dict
from cycle import reactionnetwork
import time
import billiard as multiprocessing

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

    # parser.add_argument(
    #     "-c", "--charg",
    #     default=",",
    #     dest="charg",
    #     help="The separator to change to"
    # )

    args = parser.parse_args()

    import warnings

    warnings.filterwarnings("ignore")

    kwargs, rec = read_input(args.in_file)
    if not rec:
        print("Failed to read input")
        sys.exit()

    # determine the working folder
    wd = os.path.abspath(os.path.dirname(args.in_file))
    # if not os.path.isabs(wd):
    #     wd = os.path.abspath(args.in_file)
    #     if not wd:
    #         wd = "."
    #     else:
    #         wd = os.path.join('.', wd)
    print(f" working folder is {wd}")

    # clean up old resultes
    l = os.listdir(wd)
    for i in l:
        if os.path.isdir(f'{wd}/{i}'):
            shutil.rmtree(f'{wd}/{i}')
        else:
            if i.endswith('json'):
                os.remove(f'{wd}/{i}')

    # recording starting time
    t1 = time.time()

    # prepare reactants before cycle: do opt and sp and get basic propoties for each reactant
    pool = multiprocessing.Pool(len(rec.keys()))
    results = []
    for i in list(rec.keys()):
        results.append(pool.apply_async(prepare, (wd, i, rec), kwargs))
    pool.close()
    pool.join()

    for res in results:
        rec = update_dict(rec, res.get())

    # begin cycle
    rec = reactionnetwork(rec, wd, kwargs=kwargs)

    # recording finishing time
    t2 = time.time()
    print(f"finished with {t2 - t1}")
    sys.exit()
