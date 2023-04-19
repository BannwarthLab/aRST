import os 
import shutil

import numpy as np
import pandas as pd
import alignfunction 

#align(target_name_1:str, target_name_2:str, c_1:int, c_2:int, u_1:int, u_2:int):
if __name__ == "__main__":
    tmpdir = os.getcwd()
    alignfunction.align('radical', 'ClF', 0, 0, 1, 0)

    os.chdir(os.path.join(tmpdir,"scan"))
    
    os.system('xtb product.xyz --opt --input scan.inp > xtb.out')