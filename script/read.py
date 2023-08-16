import os
import re
import time
import traceback
from collections import Counter

import pandas as pd
from distutils.util import strtobool
# from callsys import write_xyz, creat_folder, submit_orca_sp, submit_orca_opt, \
#     submit_xtb_opt, collect_energy, submit_xtb_sp
from submit_tmp import seperate, get_charge, creat_folder, \
    write_xyz, submit_xtb_opt, get_nelectron, submit_orca_sp, submit_orca_opt, submit_xtb_sp, collect_energy
from get_molbar import do_molbar


def read_input(input):
    pre_wd = os.path.dirname(input)
    # read input content
    with open(input) as f:
        content = f.readlines()
    kwargs = {}
    rec = {}
    count = 1
    for i in content:
        if "molbar" in i:
            kwargs['molbar'] = strtobool([j for j in re.split("\t| |\n", i)
                                          if j != "" and "molbar" not in j][0])

    for i in content:
        # read reactants
        if "coord" in i:
            con_l = [j for j in re.split("\t| |\n", i) if j != "" and "coord" not in j]
            rec['1/' + str(count)] = {}
            filename, charg, mul = con_l[0].split(',')
            coord = pd.read_table(f'{pre_wd}/{filename}', header=None, skiprows=2, sep='\s+')
            rec['1/' + str(count)]['coord'] = coord
            rec['1/' + str(count)]['charge'] = int(charg)
            rec['1/' + str(count)]['m'] = int(mul)
            count += 1

        # read calculation setting
        if "intercands" in i:
            kwargs['intercands'] = int([j for j in re.split("\t| |\n", i)
                                        if j != "" and "intercands" not in j][0])

        if "initial_opt" in i:
            kwargs['initial_opt'] = strtobool([j for j in re.split("\t| |\n", i)
                                               if j != "" and "initial_opt" not in j][0])
        if "dft_level" in i:
            kwargs['dft_level'] = [j for j in re.split("\t| |\n", i)
                                   if j != "" and "dft_level" not in j][0]
        if "elim" in i:
            kwargs['autoscreen'] = True
            kwargs['elim'] = int([j for j in re.split("\t| |\n", i)
                                  if j != "" and "elim" not in j][0])
        else:
            kwargs['autoscreen'] = False

        if "cycle" in i:
            kwargs['cycle'] = int([j for j in re.split("\t| |\n", i)
                                   if j != "" and "cycle" not in j][0])
        if "xtbwd" in i:
            kwargs['xtbwd'] = [j for j in re.split("\t| |\n", i)
                               if j != "" and "xtbwd" not in j][0]
        if "alpb" in i:
            kwargs['alpb'] = [j for j in re.split("\t| |\n", i)
                              if j != "" and "alpb" not in j][0]
        if "reactor_like" in i:
            kwargs['reactor_like'] = strtobool([j for j in re.split("\t| |\n", i)
                                                if j != "" and "reactor_like" not in j][0])
        if "hydrogen" in i:
            kwargs['hydrogen'] = strtobool([j for j in re.split("\t| |\n", i)
                                           if j != "" and "hydrogen" not in j][0])  # for disssociation consider H
            # element

    return kwargs, rec


def prepare(wd, reactant, rec, **kwargs):
    i = reactant
    initial_opt = kwargs.get('initial_opt', True)
    level = kwargs.get('dft_level', 'r2SCAN-3c')
    alpb = kwargs.get('alpb', None)

    # xtb opt
    creat_folder(f'{wd}/{i}/xtbopt')  # forlder xtbopt not only for opt, also sp
    coord = rec[i]['coord']
    nat = coord.shape[0]
    charg, m = rec[i]['charge'], rec[i]['m']
    write_xyz(coord, 'coord')
    if nat > 1:
        submit_xtb_opt(f'{wd}/{i}/xtbopt', charg, alpb)
        while not os.path.isfile(f'{wd}/{i}/xtbopt/xtb.out'):
            time.sleep(5)
        else:
            rec = collect_energy(rec, f'{wd}/{i}/xtbopt', f'{i}', "gfn2")
            if os.path.isfile(f'{wd}/{i}/xtbopt/xtbopt.xyz'):
                coord = pd.read_table(f'{wd}/{i}/xtbopt/xtbopt.xyz', header=None,
                                      skiprows=2, sep='\s+')
                rec[i]['coord'] = coord

    else:
        submit_xtb_sp(f'{wd}/{i}/xtbopt', charg, alpb)
        while not os.path.isfile(f'{wd}/{i}/xtbopt/xtb.out'):
            time.sleep(2)
        else:
            rec = collect_energy(rec, f'{wd}/{i}/xtbopt', f'{i}', "gfn2")

    # dft
    creat_folder(f'{wd}/{i}/orca')
    write_xyz(coord, 'coord')
    if initial_opt and nat > 1:
        if charg != 0:
            submit_orca_opt(f'{wd}/{i}/orca', level, charg, m, alpb)
        else:
            submit_orca_opt(f'{wd}/{i}/orca', level, charg, m)
    else:
        if charg != 0:
            submit_orca_sp(f'{wd}/{i}/orca', level, charg, m, alpb)
        else:
            submit_orca_sp(f'{wd}/{i}/orca', level, charg, m)

    # molbar
    rec[i]['mb'] = do_molbar(coord, wd)

    # collecting dft results
    while not os.path.isfile(f'{wd}/{i}/orca/orca.out'):
        time.sleep(5)
    else:
        if os.path.isfile(f'{wd}/{i}/orca/orca.xyz'):
            rec[i]['coord'] = pd.read_table(f'{wd}/{i}/orca/orca.xyz', header=None, skiprows=2, sep='\s+')
        else:
            rec[i]['coord'] = pd.read_table(f'{wd}/{i}/orca/coord.xyz', header=None, skiprows=2, sep='\s+')

        if os.path.isfile(f'{wd}/{i}/orca/ENERGY'):
            with open(f'{wd}/{i}/orca/ENERGY', 'r') as f:
                sp_eh = float(f.readlines()[0])
            rec[i]['sp_dft'] = sp_eh
        else:
            rec[i]['sp_dft'] = float('inf')

    return rec


def update_dict(rec_old, rec_new):
    if rec_new:
        for i in list(rec_new.keys()):
            if i in rec_old.keys() and rec_new[i]:
                if type(rec_new[i]) == list:
                    rec_old[i] = rec_new[i]
                else:
                    for j in list(rec_new[i].keys()):
                        if j not in rec_old[i].keys():
                            rec_old[i][j] = {}
                            rec_old[i][j] = rec_new[i][j]
                        else:
                            rec_old[i][j] = rec_new[i][j]
            else:
                if type(rec_new[i]) == list:
                    rec_old[i] = rec_new[i]
                else:
                    rec_old[i] = {}
                    rec_old[i] = rec_new[i]
    return rec_old


def existreaction(reaction_rec, new_reaction):  # todo Database Quick Search
    reactants_l = [reaction_rec[k][0] for k in list(reaction_rec.keys())]
    trigger = False
    for i in reactants_l:
        if Counter(i) == Counter(new_reaction):
            trigger = True
    return trigger
