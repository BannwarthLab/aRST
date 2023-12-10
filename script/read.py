import builtins
import os
import re
import shutil
import sys
import time
import traceback
from dictionaryInspector import DictionaryInspector, convert_json2dict
import toml
from toolbox import get_nelectron, create_folder, get_energy, write_xyz
from callsys import CallxTB, CallORCA
import json
import pandas as pd
from struc import StrucOperation


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
    if input_d.get('searching_setting').get('continue_sesarch', False) and os.path.isfile(f'{wd}/rec.json'):
        # continuing the last search
        with open(f'{wd}/rec.json', 'r') as f:
            content = json.load(f)
        rec = convert_json2dict(content)
    else:
        # clean up old resultes
        l = os.listdir(wd)
        for i in l:
            if os.path.isdir(f'{wd}/{i}'):
                shutil.rmtree(f'{wd}/{i}')
            elif i.endswith('out') or i.endswith('json'):
                os.remove(f'{wd}/{i}')

        # read general strucal info
        geom_d = input_d['geom']
        prefix = os.path.dirname(path)
        if not os.path.isabs(path):
            prefix = os.path.join('.', prefix)

        rec = {}
        struc_inspector = DictionaryInspector(rec)
        count = 1
        pcount = 1
        for i in geom_d.keys():
            filename = input_d['geom'][i]['name']
            if "reactant" in i:
                struc_id = f"{count:04d}"
                count += 1

            elif "product" in i:
                struc_id = f"refp{pcount:04d}"
                pcount += 1
            else:
                continue
            struc_inspector.set_value(f"{struc_id}.struc_info.general.coord",
                                      pd.read_table(f'{prefix}/{filename}', header=None, skiprows=2, sep='\s+'))
            struc_inspector.set_value(f"{struc_id}.struc_info.general.charge", input_d['geom'][i]['charge'])
            struc_inspector.set_value(f"{struc_id}.struc_info.general.multiplicity", input_d['geom'][i]['multiplicity'])
    return rec, input_d


def prepare(wd, struc_id, struc_inspector, **kwargs):
    coord = struc_inspector.get_value(f"{struc_id}.struc_info.general.coord")
    charge = struc_inspector.get_value(f"{struc_id}.struc_info.general.charge")
    m = struc_inspector.get_value(f"{struc_id}.struc_info.general.multiplicity")

    tmprec = {}
    wd = create_folder(os.path.join(wd, "struc"))
    # gfn2 opt
    xtb_setting = kwargs.get('xtb_setting')
    extracmd = xtb_setting.get('extracmd', None)
    level = xtb_setting.get('xtblevel', 'gfn2')
    alpb = xtb_setting.get('alpb', 'h2o')
    etemp = xtb_setting.get("etemp", str(300))
    struc_callxtb = CallxTB(level, "coord", charge, m - 1, alpb, etemp)
    output_folder = os.path.join(wd, f"{struc_id}/xtbopt")
    create_folder(output_folder)
    write_xyz(coord, "coord")

    if coord.shape[0] > 1:
        rc = struc_callxtb.callopt(extracmd)
        if rc == 0:
            coord = pd.read_table(os.path.join(output_folder, 'xtbopt.xyz'), header=None, skiprows=2, sep='\s+')
    else:
        rc = struc_callxtb.callsp()

    if rc == 0:
        nel = get_nelectron(os.path.join(output_folder, 'xtb.out'))
        extb = get_energy(os.path.join(output_folder, 'ENERGY'))
    else:
        extb = float("inf")
        nel = float("inf")

    # get mb
    mb = StrucOperation(coord, m, charge).get_molbarcode()
    formula = mb.split('|')[2]

    # dft
    orca_setting = kwargs.get('orca_setting')
    method = orca_setting.get("dft_level", None)
    eorca = float("inf")
    if method != "gfn2":
        pcm = orca_setting.get("pcm", None)
        output_folder = os.path.join(wd, f"{struc_id}/orca")
        create_folder(output_folder)
        write_xyz(coord, "coord")
        struc_callorca = CallORCA(method, "coord", charge, m, pcm)
        if coord.shape[0] > 1 and method is not None and orca_setting.get("initial_dftopt", True):
            rc = struc_callorca.callcalculation(caltyp="sp OPT")
            if rc == 0:
                coord = pd.read_table(os.path.join(output_folder, 'orca.xyz'), header=None, skiprows=2, sep='\s+')
        else:
            rc = struc_callorca.callcalculation()

        if rc != 0:
            print(f"{struc_id} calling orca wrong")
        else:  # Reduce storage space
            eorca = get_energy(os.path.join(output_folder, "ENERGY"))
            extensions_to_remove = ['.gbw', '.densities', '.opt']
            for filename in os.listdir('.'):
                if any(filename.endswith(ext) for ext in extensions_to_remove):
                    os.remove(filename)

    # collect data
    struc_inspector = DictionaryInspector(tmprec)
    struc_inspector.set_value(f"{struc_id}.struc_info.general.coord", coord)
    struc_inspector.set_value(f"{struc_id}.struc_info.general.nel", nel)
    struc_inspector.set_value(f"{struc_id}.struc_info.general.mb", mb)
    struc_inspector.set_value(f"{struc_id}.struc_info.general.formula", formula)
    struc_inspector.set_value(f"{struc_id}.struc_info.energy.sp_gfn2", extb)
    struc_inspector.set_value(f"{struc_id}.struc_info.energy.sp_dft", eorca)
    if "refp" in struc_id:
        struc_inspector.set_value(f"{struc_id}.reaction_path.as_product", ['ref'])
    else:
        struc_inspector.set_value(f"{struc_id}.reaction_path.as_reactant", ['1'])

    # trv,solv
    if orca_setting.get("reaction_energy", 'sp') == 'g' or orca_setting.get("reaction_energy", 'sp') == 'G':
        output_folder = os.path.join(wd, f"{struc_id}/trv")
        create_folder(output_folder)
        write_xyz(coord, "coord")
        rc = struc_callxtb.calltrv()
        etrv = get_energy(os.path.join(output_folder, 'GTRV')) if rc == 0 else None

        output_folder = os.path.join(wd, f"{struc_id}/solv")
        create_folder(output_folder)
        write_xyz(coord, "coord")
        rc = struc_callxtb.callsolv()
        esolv = get_energy(os.path.join(output_folder, 'GSOLV')) if rc == 0 else None

        struc_inspector.set_value(f"{struc_id}.struc_info.energy.sp_trv", etrv)
        struc_inspector.set_value(f"{struc_id}.struc_info.energy.sp_solv", esolv)
        if etrv is not None and esolv is not None:
            eG = etrv+esolv+extb if method == "gfn2" else etrv+esolv+eorca
            struc_inspector.set_value(f"{struc_id}.struc_info.energy.G", eG)
    return tmprec
