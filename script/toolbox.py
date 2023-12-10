import os
import re
import shutil
import string
import subprocess
import json
import copy
import random
import traceback
import math
import numpy as np
import pandas as pd


def create_folder(wd):
    wd = os.path.abspath(wd)
    if not os.path.exists(wd) and not os.path.isdir(wd):
        try:
            os.makedirs(wd)
        except:
            pass
    os.chdir(wd)
    return wd

def rm_folderORfile(wd):
    wd = os.path.abspath(wd)
    if os.path.isdir(wd):
        try:
            shutil.rmtree(f'{wd}', ignore_errors=True)
            os.rmdir(f'{wd}')
        except:
            pass
    elif os.path.isfile(wd):
        try:
            os.remove(wd)
        except:
            pass




def creat_randomfolder(wd):
    random_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
    if not os.path.isdir(f'{wd}/{random_name}'):
        os.makedirs(f'{wd}/{random_name}')
    else:
        try:
            shutil.rmtree(f'{wd}/{random_name}', ignore_errors=True)
            os.rmdir(f'{wd}/{random_name}')
        except:
            pass
        os.makedirs(f'{wd}/{random_name}')
    os.chdir(f'{wd}/{random_name}')
    return f'{wd}/{random_name}'


def write_xyz(xyz, name):
    pwd = os.getcwd()
    # xyz.round(10)
    nat = xyz.shape[0]
    with open(f'{pwd}/{name}.xyz', 'w') as f:
        f.write(f'{nat}\n')
        f.write('\n')
        for i in range(nat):
            at, x, y, z = xyz.iloc[i, 0], xyz.iloc[i, 1], xyz.iloc[i, 2], xyz.iloc[i, 3]
            f.write('%2s %20.12f %20.12f %20.12f\n' % (at, x, y, z))
        f.write('\n')


def seperate():
    command = f'excprep coord.xyz > e.out'
    subprocess.run(command, shell=True)
    return


def get_charge(coord_new, coord_old, chargefile):
    context = pd.read_table(chargefile, header=None, sep='\s+')

    nat_l = []
    coord_new = coord_new.round(4)
    coord_old = coord_old.round(4)
    try:
        for i in range(coord_new.shape[0]):
            new = [coord_new.iloc[i, 1], coord_new.iloc[i, 2], coord_new.iloc[i, 3]]
            for j in range(coord_old.shape[0]):
                old = [coord_old.iloc[j, 1], coord_old.iloc[j, 2], coord_old.iloc[j, 3]]
                if new == old:
                    nat_l.append(j)
    except IndexError:
        print(chargefile)
    charge = 0
    for i in nat_l:
        charge += float(context.iloc[i, 0])

    # tmp set homogeneous dissociation
    a = int(math.modf(charge)[1])
    if abs(charge - a) >= 0.5:
        return int(np.round(charge))
    else:
        return a


def get_uhf_from_excprep(excprepoutwd, fnr):
    uhf = {}
    with open(f'{excprepoutwd}', 'r+') as f:
        lines = f.readlines()
    for i in range(1, fnr + 1):
        for j in lines:
            if f"{str(i)}:" in j:
                if "UHF" not in j:
                    uhf[i] = 0
                else:
                    uhf[i] = int(re.findall(r'\-?\d+\.?\d*', str(re.findall(re.compile(r'[(](.*?)[)]', re.S), j)))[0])
    return uhf


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


def get_energy(ENERGYFILE):
    speh = float("inf")
    if os.path.isfile(f'{ENERGYFILE}'):
        with open(f'{ENERGYFILE}', 'r') as f:
            speh = float(f.readlines()[0])
    return speh



def boltzmanp(a):
    # a: gap or wb(wboe)

    ev2j = 1.602176634 * pow(10, -19)
    return np.exp(-a * ev2j / (1.38 * pow(10, -23) * 1000))


def lowestgeom(infile):
    # input is .log file
    with open(infile, 'r') as f:
        f.seek(0, 0)
        first = f.readline()
        atnr = int(first.split('\n')[0])
        f.seek(0, 0)
        d = {}
        test = f.readlines()
        for i in range(0, len(test), atnr + 2):
            tmp = [test[j] for j in range(i + 1, i + atnr + 2)]
            try:
                energy = float(re.findall(r'\-?\d+\.?\d*', tmp[0].split(':')[1])[0])
            except:
                energy = 'inf'
            d[energy] = tmp
        a = min(d, key=d.get)
        outlist = [str(atnr) + '\n'] + d[a]
        return outlist


def sum_energy(fragments,level,struc_inspector):
    energy = None
    tmp = [struc_inspector.get_value(f"{k}.struc_info.energy.{level}") for k in fragments]
    print("fragments:",fragments)
    print("energy:",tmp)
    if not any(element is None for element in tmp):
        energy= sum(tmp)
    return energy

def get_reactionenergy(reactants, products, level,struc_inspector):
    ereaction = None
    ereactants = sum_energy(reactants, level,struc_inspector)
    eproducts = sum_energy(products, level,struc_inspector)
    if ereactants is not None and eproducts is not None:
        ereaction = eproducts - ereactants
    return ereaction

# def collect_sp_xyz(rec, reckey, wd, typ):
#     speh = float('inf')
#     if typ == 'gfn2':
#         if os.path.isfile(f'{wd}/ENERGY'):
#             with open(f'{wd}/ENERGY', 'r') as f:
#                 speh = float(f.readlines()[0])
#         if 'sp_gfn2' not in rec[reckey].keys() or rec[reckey]['sp_gfn2'] > speh:
#             rec[reckey]['sp_gfn2'] = speh
#             if os.path.isfile(f'{wd}/xtbopt.xyz'):
#                 rec[reckey]['coord'] = pd.read_table(f'{wd}/xtbopt.xyz', header=None, skiprows=2, sep='\s+')
#             else:
#                 rec[reckey]['coord'] = pd.read_table(f'{wd}/coord.xyz', header=None, skiprows=2, sep='\s+')
#             nat = rec[reckey]['coord'].shape[0]
#             if nat == 1 and rec[reckey]['coord'].iloc[0, 0] in ["H", "h"]:
#                 rec[reckey]['sp_gfn2'] = 0
#     else:
#         if os.path.isfile(f'{wd}/ENERGY'):
#             with open(f'{wd}/ENERGY', 'r') as f:
#                 speh = float(f.readlines()[0])
#         if 'sp_dft' not in rec[reckey].keys() or rec[reckey]['sp_dft'] > speh:
#             rec[reckey]['sp_dft'] = speh
#             if os.path.isfile(f'{wd}/orca.xyz'):
#                 rec[reckey]['coord'] = pd.read_table(f'{wd}/orca.xyz', header=None, skiprows=2, sep='\s+')
#             else:
#                 rec[reckey]['coord'] = pd.read_table(f'{wd}/coord.xyz', header=None, skiprows=2, sep='\s+')
#             nat = rec[reckey]['coord'].shape[0]
#             if nat == 1 and rec[reckey]['coord'].iloc[0, 0] in ["H", "h"] and rec[reckey]['charge'] == 1:
#                 rec[reckey]['sp_dft'] = 0
#                 rec[reckey]['m'] = 0
#             else:
#                 with open(f'{wd}/orca.inp', 'r') as f:
#                     string = f.read()
#                 pattern = r'\* xyzfile -?\d+ (\d+) coord\.xyz'
#                 m = re.search(pattern, string).group(1)
#                 rec[reckey]['m'] = int(m)
#     return rec


# def collect_g_xyz(rec, reckey, wd):
#     speh1 = 0
#     speh2 = 0
#     if os.path.isfile(f'{wd}/gfn2_trv/GTRV'):
#         with open(f'{wd}/gfn2_trv/GTRV', 'r') as f:
#             speh1 = float(f.readlines()[0])
#     if os.path.isfile(f'{wd}/gfn2_solv/GSOLV'):
#         with open(f'{wd}/gfn2_solv/GSOLV', 'r') as f:
#             speh2 = float(f.readlines()[0])
#     sp = rec[reckey]['sp_dft']
#     rec[reckey]['g'] = sp + speh1 + speh2
#     rec[reckey]['trv'] = speh1
#     rec[reckey]['solv'] = speh2
#     return rec
