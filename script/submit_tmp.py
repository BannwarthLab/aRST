import os
import re
import shutil
import subprocess
import json
import copy
import random
import time

import numpy as np
import pandas as pd
from tempfile import mkdtemp


def remove_file(old_path, new_path):
    filelist = os.listdir(old_path)

    for file in filelist:
        src = os.path.join(old_path, file)
        dst = os.path.join(new_path, file)

        shutil.move(src, dst)


def submit_xtb_sp(wd, charg, alpb=None):
    if alpb:
        command = f'xtb coord.xyz -c {charg}  --alpb {alpb} --etemp 1000 > xtb.out ' \
                  f'&& tac xtb.out | grep -m 1 \'TOTAL ENERGY\' |awk \'{{print $4}}\' > ENERGY \n'
    else:
        command = f'xtb coord.xyz -c {charg} --etemp 1000  > xtb.out && ' \
                  f'tac xtb.out | grep -m 1 \'TOTAL ENERGY\' |awk \'{{print $4}}\' > ENERGY \n'
    tmpd = mkdtemp()
    remove_file(wd, tmpd)
    os.chdir(tmpd)
    subprocess.run(command, shell=True)
    while not os.path.isfile(f"{tmpd}/xtb.out"):
        time.sleep(1)
    else:
        remove_file(tmpd, wd)
        os.chdir(wd)
        if os.path.isdir(tmpd):
            shutil.rmtree(tmpd)


def submit_xtb_opt(wd, charg, alpb=None):
    if alpb:
        command = f'xtb coord.xyz --opt -c {charg}  --alpb {alpb} --etemp 1000 > xtb.out ' \
                  f'&& tac xtb.out | grep -a \'TOTAL ENERGY\' |awk \'{{print $4}}\' > ENERGY \n'
    else:
        command = f'xtb coord.xyz --opt -c {charg} --etemp 1000  > xtb.out && ' \
                  f'tac xtb.out | grep -a \'TOTAL ENERGY\' |awk \'{{print $4}}\' > ENERGY \n'
    tmpd = mkdtemp()
    remove_file(wd, tmpd)
    os.chdir(tmpd)
    subprocess.run(command, shell=True)
    while not os.path.isfile(f"{tmpd}/xtb.out"):
        time.sleep(1)
    else:
        remove_file(tmpd, wd)
        os.chdir(wd)
        if os.path.isdir(tmpd):
            shutil.rmtree(tmpd)


def submit_scan(wd, charg, alpb=None, **kw):
    concerted_trigger = kw.get('concerted_trigger', False)
    kw = kw.get('kw')
    at1 = kw['at1']
    at2 = kw['at2']
    d1 = kw['d1']
    d2 = kw['d2']
    if concerted_trigger:
        at3 = kw['at3']
        at4 = kw['at4']
        d_34_start = kw['d_34_start']
        d_34_end = kw['d_34_end']
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant=1\n'
                     f'$scan\n'
                     'mode=concerted\n'
                     f'distance: {at3},{at4},{d_34_start};{d_34_start},{d_34_end},30\n'
                     f'distance: {at1},{at2},{d1};{d1},{d2},30\n$end')
    else:
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant=1\n'
                     f'$scan\n'
                     f'distance: {at1},{at2},{d1};{d1},{d2},30\n$end')

    if alpb:
        command = f'xtb coord.xyz --opt  -c {charg}  --alpb {alpb} --etemp 1000 --input scan.inp > xtb.out\n'
    else:
        command = f'xtb coord.xyz --opt  -c {charg} --etemp 1000 --input scan.inp > xtb.out\n'
    tmpd = mkdtemp()
    remove_file(wd, tmpd)
    os.chdir(tmpd)
    print(command)
    subprocess.run(command, shell=True)
    while not os.path.isfile(f"{tmpd}/xtb.out"):
        time.sleep(1)
    else:
        remove_file(tmpd, wd)
        os.chdir(wd)
        if os.path.isdir(tmpd):
            shutil.rmtree(tmpd)


# f" %pal nprocs 8 end \n" \
def submit_orca(wd, method, charg, m, nat, alpb=None):
    if nat > 1:
        caltyp = 'sp OPT'
    else:
        caltyp = 'sp'
    solmode = ""
    if charg:
        if alpb is None:
            print(f"please use CPCM mode for {wd}")
        else:
            if alpb == "h2o" or "H2O":
                alpb = "water"
            solmode = f"CPCM({alpb})"
    if m == 1:
        caltyp2 = "RKS"
    else:
        caltyp2 = "UKS"
    cmd = f"! {method} {caltyp} {caltyp2} {solmode}\n" \
          f"* xyzfile {charg} {m} coord.xyz "
    with open(f'{wd}/orca.inp', 'w') as f:
        f.write(cmd)

    ######## if running in debug hpz107
    tmpd = mkdtemp()
    remove_file(wd, tmpd)
    os.chdir(tmpd)

    cmd = f'$ORCA orca.inp  > orca.out && tac orca.out | grep -m 1 \'FINAL SINGLE\' |awk \'{{print $5}}\' > ENERGY'
    subprocess.run(cmd, shell=True)
    while not os.path.isfile(f"{tmpd}/orca.out"):  # double check if shell returns before finish
        time.sleep(1)
    else:
        remove_file(tmpd, wd)
        os.chdir(wd)
        if os.path.isdir(tmpd):
            shutil.rmtree(tmpd)





def creat_folder(wd):
    if not os.path.isdir(f'{wd}'):
        os.makedirs(f'{wd}')
    os.chdir(f'{wd}')
    return


def creat_randomfolder(wd):
    random_name = "tmp_" + str(random.randint(1, 999))
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
    charg = int(np.round(charg, 0))

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
            rec[f'{name}'][f'sp_{type}'] = sp_eh
        except:
            print(f"{name}, sp_{type} energy is wrong")
            rec[f'{name}'][f'sp_{type}'] = 0
    else:
        try:
            nat = rec[f'{name}']['coord'].shape[0]
            if nat == 1 and rec[f'{name}']['coord'].iloc[0, 0] in ["H", "h"]:
                rec[f'{name}'][f'sp_{type}'] = 0
            else:
                rec[f'{name}'][f'sp_{type}'] = float('inf')
        except:
            rec[f'{name}'][f'sp_{type}'] = float('inf')
    return rec


def write_rec(rec, td, react_l=None):
    crec = copy.deepcopy(rec)
    if react_l is not None:
        react_l = [react_l][0]
        tmp = []
        tmp_l = []
        for i in list(crec.keys()):
            if 'mb' in crec[i].keys():
                mb = crec[i]['mb']
                if mb not in tmp:
                    tmp.append(mb)
                    tmp_l.append(i)

        for i in list(crec.keys()):
            if i not in tmp_l:
                del crec[i]

    for name in crec:
        for i in list(crec[name].keys()):
            if isinstance(crec[name][i], pd.DataFrame):
                crec[name][i] = crec[name][i].to_json()

    if react_l:
        # clean up old rec
        if os.path.isfile(f'{td}/reactor.json'):
            os.remove(f'{td}/reactor.json')
        with open(f'{td}/reactor.json', "a+") as f:
            f.write(json.dumps(crec))
            f.write('\n')
    else:
        # clean up old rec
        if os.path.isfile(f'{td}/rec.json'):
            os.remove(f'{td}/rec.json')
        with open(f'{td}/rec.json', "a+") as f:
            f.write(json.dumps(crec))
            f.write('\n')


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
