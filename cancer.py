import os 
import shutil

import numpy as np
import pandas as pd
import RS_singlestep

#align(target_name_1:str, target_name_2:str, c_1:int, c_2:int, u_1:int, u_2:int):

if __name__ == "__main__":
    reaction_type = 1
    tmpdir = os.getcwd()
    if reaction_type == 0:
        RS_singlestep.intermolecular_reaction('radical', 'ClF', 0, 0, 1, 0)
        os.chdir(os.path.join(tmpdir,"scan"))
        os.system('xtb product.xyz --opt --input scan.inp > xtb.out')
    elif reaction_type == 1:
        RS_singlestep.call_dissociation('methane', 0, 0)
        os.chdir(os.path.join(tmpdir,"scan"))
        os.system('xtb methane.xyz --opt --input scan.inp > xtb.out')
