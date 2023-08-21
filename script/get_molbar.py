import shutil
import traceback
import time
import pandas as pd
from molbar.barcode import get_molbar
import os
from operation import write_xyz, creat_randomfolder


def do_molbar(coord, wd):
    out = []
    if isinstance(coord, pd.DataFrame):
        pwd = os.getcwd()
        fname = creat_randomfolder(wd)
        write_xyz(coord, 'coord')
        # out = MolBar().calculate([f'{wd}/tmp/coord.xyz'], threads=8)[0][1] # molbar 0.6.1
        out = get_molbar([f'{wd}/{fname}/coord.xyz'], threads=8)[0][1]
        while not out:
            time.sleep(0.5)
        else:
            os.chdir(pwd)
            shutil.rmtree(f'{wd}/{fname}', ignore_errors=True)
    return out
