import builtins
import os
import re
import sys
import time
import pandas as pd
from collections import Counter
from get_molbar import do_molbar
from callsys import call_orca_dft, call_xtb
from operation import creat_folder, write_xyz, collect_energy


def Convert(s):
    if s == "True" or s == "TRUE":
        return True
    else:
        return False


def read_input(path):
    # read input content
    with open(path) as f:
        content = f.readlines()
    kwargs = {}
    rec = {}
    count = 1
    pcount = 1
    for i in content:
        if "molbar" in i:
            kwargs['molbar'] = Convert([j for j in re.split("\t| |\n", i)
                                        if j != "" and "molbar" not in j][0])

    prefix = os.path.dirname(path)
    if not os.path.isabs(path):
        prefix = os.path.join('.', prefix)

    for i in content:
        # read reactants
        if "coord" in i:
            if "pcoord" in i:  # for test
                con_l = [j for j in re.split("\t| |\n", i) if j != "" and "pcoord" not in j]
                filename, charg, mul = con_l[0].split(',')
                pcoord = pd.read_table(f'{prefix}/{filename}', header=None, skiprows=2, sep='\s+')
                rec['refproduct/' + str(pcount)] = {}
                rec['refproduct/' + str(pcount)]['coord'] = pcoord
                rec['refproduct/' + str(pcount)]['charge'] = int(charg)
                rec['refproduct/' + str(pcount)]['m'] = int(mul)
                if 'refproduct/' + str(pcount) not in rec.keys():
                    print("careful! ref products didn't recorded")
                pcount += 1
            else:
                con_l = [j for j in re.split("\t| |\n", i) if j != "" and "coord" not in j]
                rec['1/' + str(count)] = {}
                filename, charg, mul = con_l[0].split(',')
                coord = pd.read_table(f'{prefix}/{filename}', header=None, skiprows=2, sep='\s+')
                rec['1/' + str(count)]['coord'] = coord
                rec['1/' + str(count)]['charge'] = int(charg)
                rec['1/' + str(count)]['m'] = int(mul)
                count += 1

        # read calculation setting
        if "intercands" in i:
            kwargs['intercands'] = int([j for j in re.split("\t| |\n", i)
                                        if j != "" and "intercands" not in j][0])
        if "discands" in i:
            kwargs['discands'] = int([j for j in re.split("\t| |\n", i)
                                      if j != "" and "discands" not in j][0])

        if "initial_opt" in i:
            kwargs['initial_opt'] = Convert([j for j in re.split("\t| |\n", i)
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
            kwargs['reactor_like'] = Convert([j for j in re.split("\t| |\n", i)
                                              if j != "" and "reactor_like" not in j][0])
        if "hydrogen" in i:
            kwargs['hydrogen'] = Convert([j for j in re.split("\t| |\n", i)
                                          if j != "" and "hydrogen" not in j][0])  # for disssociation consider H
            # element
        if "platform" in i:
            kwargs['platform'] = int([j for j in re.split("\t| |\n", i)
                                      if j != "" and "platform" not in j][0])
            # mode 1 : only run in debbug node: call program with commandline, but doing cal in /tmp
            # mode 2 : run py script locally and submit each process with slurm
            # NONE/0: default mode, sbatch this py scrit
        if "constrain_strength" in i:
            kwargs['constrain_strength'] = int([j for j in re.split("\t| |\n", i)
                                                if j != "" and "constrain_strength" not in j][0])
            # constrain_strength in xtbscan
    return kwargs, rec


def prepare(wd, i, rec, **kwargs):
    initial_opt = kwargs.get('initial_opt', True)
    level = kwargs.get('dft_level', None)
    alpb = kwargs.get('alpb', None)
    mode = kwargs.get('platform', 0)
    if not level:
        typ = 'gfn2'
    else:
        typ = 'dft'

    # xtb opt
    creat_folder(f'{wd}/{i}/xtbopt')  # forlder xtbopt not only for opt, also sp
    coord = rec[i]['coord']
    nat = coord.shape[0]
    charg, m = rec[i]['charge'], rec[i]['m']
    write_xyz(coord, 'coord')
    call_xtb(f'{wd}/{i}/xtbopt', charg, nat, uhf=m - 1, alpb=alpb, mode=mode)
    if os.path.isfile(f'{wd}/{i}/xtbopt/xtb.out'):
        rec = collect_energy(rec, f'{wd}/{i}/xtbopt', f'{i}', "gfn2")
        if os.path.isfile(f'{wd}/{i}/xtbopt/xtbopt.xyz'):
            coord = pd.read_table(f'{wd}/{i}/xtbopt/xtbopt.xyz', header=None,
                                  skiprows=2, sep='\s+')
            rec[i]['coord'] = coord
    else:
        print("call xtb wrong")
        sys.exit()

    # molbar
    rec[i]['mb'] = do_molbar(coord, wd)

    # dft
    if typ == 'dft':
        creat_folder(f'{wd}/{i}/orca')
        write_xyz(coord, 'coord')
        if not initial_opt:
            call_orca_dft(f'{wd}/{i}/orca', level, charg, m, nat, opt=False, alpb=alpb, mode=mode)
        else:
            call_orca_dft(f'{wd}/{i}/orca', level, charg, m, nat, alpb=alpb, mode=mode)
        if os.path.isfile(f'{wd}/{i}/orca/orca.out'):
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
        else:
            print("call orca wrong")
    return rec


def update_dict(rec_old, rec_new):
    if rec_new:
        for i in list(rec_new.keys()):
            if type(rec_new[i]) == list:
                rec_old[i] = rec_new[i]
            else:
                if i in rec_old.keys():
                    for j in list(rec_new[i].keys()):
                        if j not in rec_old[i].keys():
                            rec_old[i][j] = {}
                            rec_old[i][j] = rec_new[i][j]
                        else:
                            rec_old[i][j] = rec_new[i][j]
                else:
                    rec_old[i] = {}
                    rec_old[i] = rec_new[i]
    return rec_old


def existreaction(reaction_rec, new_reaction, typ):  # todo Database Quick Search
    if typ == "r":
        reactants_l = [reaction_rec[k][0] for k in list(reaction_rec.keys())]
    else:
        reactants_l = [reaction_rec[k][1] for k in list(reaction_rec.keys())]
        reactants_l.remove(new_reaction)
    trigger = False
    for i in reactants_l:
        if Counter(i) == Counter(new_reaction):
            trigger = True
    return trigger
