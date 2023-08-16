import os
import traceback
import numpy as np
import pandas as pd
import ast
import random
import math

hydrogen = True
wboe = pd.read_csv(f'/Users/cy/jtc/Meow/test4/1/1/ml/wboe', header=None, sep='\s+')
if os.path.getsize(f'/Users/cy/jtc/Meow/test4/1/1/ml/wboe') != 0:
    mlfeatures = pd.read_csv(f'/Users/cy/jtc/Meow/test4/1/1/ml/ml_feature.csv')
    tmp = {}
    for i in range(wboe.shape[0]):
        dagp = wboe.iloc[i, 3]
        at_1, at_2 = wboe.iloc[i, 0], wboe.iloc[i, 1]
        ele_1 = mlfeatures.loc[at_1 - 1, 'Atom']
        ele_2 = mlfeatures.loc[at_2 - 1, 'Atom']
        if not hydrogen and (int(ele_1) == 1 or int(ele_2) == 1):
            pass
        else:
            tmp[dagp] = [ at_1, at_2]
    if not tmp:  # only h-bonds in struc
        for i in range(wboe.shape[0]):
            dagp = wboe.iloc[i, 3]
            at_1, at_2 = wboe.iloc[i, 0], wboe.iloc[i, 1]

            tmp[dagp] = [at_1, at_2]
    print(tmp)
