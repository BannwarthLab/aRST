import os
import re
import shutil
import subprocess
import json
import copy
import random
import traceback

import numpy as np
import pandas as pd


def creat_folder(wd):
    if not os.path.isdir(f'{wd}'):
        os.makedirs(f'{wd}')
    os.chdir(f'{wd}')
    return


def creat_randomfolder(wd):
    random_name = "tmp_" + str(random.randint(1, 99))
    if not os.path.isdir(f'{wd}/{random_name}'):
        os.makedirs(f'{wd}/{random_name}')
    else:
        shutil.rmtree(f'{wd}/{random_name}', ignore_errors=True)
        os.makedirs(f'{wd}/{random_name}')
    os.chdir(f'{wd}/{random_name}')
    return random_name


def write_xyz(xyz, name):
    pwd = os.getcwd()
    xyz.round(5)
    nat = xyz.shape[0]
    with open(f'{pwd}/{name}.xyz', 'w') as f:
        f.write(f'{nat}\n')
        f.write('\n')
        for i in range(nat):
            at, x, y, z = xyz.iloc[i, 0], xyz.iloc[i, 1], xyz.iloc[i, 2], xyz.iloc[i, 3]
            f.write('%2s %10.5f %10.5f %10.5f\n' % (at, x, y, z))
        f.write('\n')


def xtb_ml(xtbwd, charg, alpb=None):
    if alpb is not None:
        alpb = [alpb][0]
        cmd = f'{xtbwd} coord.xyz --ml_feature -c {charg} --alpb {alpb}  > xtb.out'
    else:
        cmd = f'{xtbwd} coord.xyz --ml_feature -c {charg}  > xtb.out'
    print(cmd)
    subprocess.run(cmd, shell=True)


def seperate():
    command = f'excprep coord.xyz > e.out'
    subprocess.run(command, shell=True)
    return


def get_charge(coord_new, coord_old, chargfile):
    nat_l = []
    coord_new = coord_new.round(4)
    coord_old = coord_old.round(4)
    for i in range(coord_new.shape[0]):
        new = [coord_new.iloc[i, 1], coord_new.iloc[i, 2], coord_new.iloc[i, 3]]
        for j in range(coord_old.shape[0]):
            old = [coord_old.iloc[j, 1], coord_old.iloc[j, 2], coord_old.iloc[j, 3]]
            if new == old:
                nat_l.append(j)
    charg = 0
    for i in nat_l:
        charg += float(chargfile.iloc[i, 0])
    charg = int(np.round(charg))

    return charg


def get_nelectron(xtboutwd):
    with open(f'{xtboutwd}', 'r') as f:
        l = []
        for line in f:
            if "# electrons" in line:
                l.append(line)
        if l:
            e = l[0]
            ne = int(re.findall(r'\-?\d+\.?\d*', e)[0])
        else:
            ne = 0
    return ne


def collect_energy(rec, wd, name, type):
    if os.path.isfile(f'{wd}/ENERGY'):
        try:
            with open(f'{wd}/ENERGY', 'r') as f:
                sp_eh = float(f.readlines()[0])
            if f'sp_{type}' in rec[f'{name}'].keys() and '_' in wd.split('/')[-1]:  # only for m= 1 or 3
                old_energy = rec[f'{name}'][f'sp_{type}']
                if sp_eh < old_energy:
                    rec[f'{name}'][f'sp_{type}'] = sp_eh
                    rec[f'{name}']['m'] = wd.split('/')[-1].split('_')[-1]
                    if os.path.isfile(f'{wd}/orca.xyz'):
                        rec[f'{name}']['coord'] = pd.read_table(f'{wd}/orca.xyz', header=None, skiprows=2, sep='\s+')
                    else:
                        rec[f'{name}']['coord'] = pd.read_table(f'{wd}/coord.xyz', header=None, skiprows=2, sep='\s+')
            else:
                rec[f'{name}'][f'sp_{type}'] = sp_eh
                if isinstance(rec[f'{name}']['m'], list) and '_' in wd.split('/')[-1]:
                    rec[f'{name}']['m'] = wd.split('/')[-1].split('_')[-1]
                    if os.path.isfile(f'{wd}/orca.xyz'):
                        rec[f'{name}']['coord'] = pd.read_table(f'{wd}/orca.xyz', header=None, skiprows=2, sep='\s+')
                    else:
                        rec[f'{name}']['coord'] = pd.read_table(f'{wd}/coord.xyz', header=None, skiprows=2, sep='\s+')
        except:
            traceback.print_exc()
            print(f"{name}, sp_{type} energy is wrong")
            rec[f'{name}'][f'sp_{type}'] = 0
    else:
        nat = rec[f'{name}']['coord'].shape[0]
        if nat == 1 and rec[f'{name}']['coord'].iloc[0, 0] in ["H", "h"]:
            rec[f'{name}'][f'sp_{type}'] = 0
        else:
            rec[f'{name}'][f'sp_{type}'] = float('inf')
    return rec


def write_rec(rec, td, filename):
    crec = copy.deepcopy(rec)
    for name in crec:
        if not type(crec[name]) == list:
            for i in list(crec[name].keys()):
                if isinstance(crec[name][i], pd.DataFrame):
                    crec[name][i] = crec[name][i].to_json()
    # clean up old rec
    if os.path.isfile(f'{td}/{filename}.json'):
        os.remove(f'{td}/{filename}.json')
    with open(f'{td}/{filename}.json', "a+") as f:
        f.write(json.dumps(crec))
        f.write('\n')
    return


def testcrop(wd, coord):
    pwd = os.getcwd()
    f = creat_randomfolder(wd)
    write_xyz(coord, 'coord')
    seperate()
    l = [k for k in os.listdir(f'{wd}/{f}') if os.path.isdir(f'{wd}/{f}/{k}') and str.isdigit(k)]
    os.chdir(pwd)
    shutil.rmtree(f'{wd}/{f}', ignore_errors=True)
    if len(l) > 1:
        return False
    else:
        return True
