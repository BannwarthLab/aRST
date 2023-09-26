import builtins
import os
import re
import shutil
import sys
import time
import traceback

import toml
from get_molbar import do_molbar
from operation import collect_sp_xyz

import json
import pandas as pd
from struc import Struc


def read_input(path):
    with open(path, 'r') as f:
        input_d = toml.load(f)
    mode = input_d.get('searching_setting').get('mode', None)
    if mode is None:
        input_d['xtb_setting']['mode'] = 0
        input_d['orca_setting']['mode'] = 0
    else:
        input_d['xtb_setting']['mode'] = mode
        input_d['orca_setting']['mode'] = mode
    homogene_dissoziation = input_d.get('xtb_setting').get('homogene', False)
    if homogene_dissoziation:
        input_d['searching_setting']['homogene'] = True

    # print(input_d)
    wd = os.path.abspath(os.path.dirname(path))
    if input_d.get('searching_setting').get('continue_sesarch', False) == True and os.path.isfile(f'{wd}/rec.json'):
        # continuing the last search
        with open(f'{wd}/rec.json', 'r') as f:
            content = json.load(f)
        rec = {}
        for struc in content.keys():
            rec[struc] = {}
            for item in list(content[struc].keys()):
                try:
                    df = pd.read_json(content[struc][item])
                    rec[struc][item] = df
                except:
                    rec[struc][item] = content[struc][item]
    else:
        # clean up old resultes
        l = os.listdir(wd)
        for i in l:
            if os.path.isdir(f'{wd}/{i}'):
                shutil.rmtree(f'{wd}/{i}')
            elif i.endswith('out') or i.endswith('json'):
                os.remove(f'{wd}/{i}')

        # read Geom
        geom_d = input_d['geom']
        prefix = os.path.dirname(path)
        if not os.path.isabs(path):
            prefix = os.path.join('.', prefix)

        rec = {}
        count = 1
        pcount = 1
        for i in geom_d.keys():
            if "reactant" in i:
                rec['1/' + str(count)] = {}
                filename = input_d['geom'][i]['name']
                charge = input_d['geom'][i]['charge']
                m = input_d['geom'][i]['multiplicity']
                coord = pd.read_table(f'{prefix}/{filename}', header=None, skiprows=2, sep='\s+')
                rec['1/' + str(count)]['coord'] = coord
                rec['1/' + str(count)]['charge'] = int(charge)
                rec['1/' + str(count)]['m'] = int(m)
                count += 1
            elif "product" in i:
                rec['refproduct/' + str(pcount)] = {}
                filename = input_d['geom'][i]['name']
                charge = input_d['geom'][i]['charge']
                m = input_d['geom'][i]['multiplicity']
                coord = pd.read_table(f'{prefix}/{filename}', header=None, skiprows=2, sep='\s+')
                rec['refproduct/' + str(pcount)]['coord'] = coord
                rec['refproduct/' + str(pcount)]['charge'] = int(charge)
                rec['refproduct/' + str(pcount)]['m'] = int(m)
                pcount += 1

    return rec, input_d


def prepare(wd, i, rec, **kwargs):
    coord = rec[i]['coord']
    charge, m = rec[i]['charge'], rec[i]['m']
    if "sp_gfn2" not in rec[i].keys():
        # xtb opt
        xtb_setting = kwargs.get('xtb_setting')
        rc = Struc(coord, m, charge).runxtb(f'{wd}/{i}/xtbopt', **xtb_setting)
        if rc != 0:
            sys.exit(f'{wd}/{i}/ xtb error')


        # collect gfn2 data
        rec = collect_sp_xyz(rec, i, f'{wd}/{i}/xtbopt/', 'gfn2')

    # molbar
    if "mb" not in rec[i].keys():
        rec[i]['mb'], rec[i]['name'] = do_molbar(coord, wd)

    # dft, sp required, opt optional
    if "sp_dft" not in rec[i].keys():
        # gfn2 overwrite option
        orca_setting = kwargs.get('orca_setting')
        if orca_setting.get('dft_level', 'r2scan-3c') != 'gfn2':
            if orca_setting.get('dft_opt', False) and orca_setting.get('initial_dftopt', True):
                # CALL ORCA OPT
                Struc(coord, m, charge).runorcasp(f'{wd}/{i}/orca', caltyp='sp OPT', **orca_setting)
            else:
                # CALL ORCA SP
                Struc(coord, m, charge).runorcasp(f'{wd}/{i}/orca', caltyp='sp', **orca_setting)

            # collect dft results
            rec = collect_sp_xyz(rec, i, f'{wd}/{i}/orca', 'dft')
    return rec


def update_dict(rec_old, rec_new):
    # update rec, reaction_rec
    if rec_new:
        for i in list(rec_new.keys()):
            try:
                if type(rec_new[i]) == list:
                    # override reaction_rec
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
            except:
                traceback.print_exc()
                pass
    return rec_old
