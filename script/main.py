import os
import shutil
import sys
from argparse import ArgumentParser
from read import read_input, prepare, update_dict
from cycle import reactionnetwork
import time
import billiard as multiprocessing

if __name__ == "__main__":
    # parser = ArgumentParser(
    #     # name of program
    #     prog="MEOW",
    #     # description
    #     description="CHECK CHECK.",
    #     # bottom text
    #     epilog="By Ying"
    # )
    # parser.parse_args()
    #
    # # add arguments
    #
    # parser.add_argument(
    #     "in_file",
    #     help="The input coord file"
    # )
    #
    #
    # parser.add_argument(
    #     "-c", "--charg",
    #     default=",",
    #     dest="charg",
    #     help="The separator to change to"
    # )
    #
    # args = parser.parse_args()
    import warnings

    warnings.filterwarnings("ignore")

    # wd = "/home/guests/ychen/meow/test3"
    wd = os.getcwd()

    ## clean up
    l = os.listdir(wd)
    for i in l:
        if os.path.isdir(f'{wd}/{i}'):
            shutil.rmtree(f'{wd}/{i}')
        else:
            if i.endswith('json'):
                os.remove(f'{wd}/{i}')

    # recoring running time
    t1 = time.time()

    # read input setting
    kwargs, rec = read_input(f'{wd}/in.inp')
    if not rec:
        print("Failed to read input")
        sys.exit()

    # do opt and sp and get basic propoties for each reactant
    pool = multiprocessing.Pool(len(rec.keys()))
    results = []

    for i in list(rec.keys()):
        results.append(pool.apply_async(prepare, (wd, i, rec), kwargs))
    pool.close()
    pool.join()

    for res in results:
        rec = update_dict(rec, res.get())

    rec = reactionnetwork(rec, wd, kwargs=kwargs)
    t2 = time.time()
    print(f"finished with {t2 - t1}")
