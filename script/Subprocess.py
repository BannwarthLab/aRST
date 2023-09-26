import builtins
import os
import shutil
import sys
import time
import traceback
import copy
import pandas as pd
from get_molbar import do_molbar
from operation import seperate, get_charge, creat_folder, write_xyz, \
    get_nelectron, testcrop, lowestgeom, get_uhf_from_excprep, collect_sp_xyz
from collections import Counter
from struc import Struc


def print(text):
    builtins.print(text)
    os.fsync(sys.stdout)


def check_xtb_scan(wd, path_i):
    del_l = []
    if not os.path.exists(f'{wd}/{path_i}/scan/xtb.out'):
        print(f'{wd}/{path_i}/scan/ failed')
        shutil.rmtree(f'{wd}/{path_i}', ignore_errors=True)
        del_l.append(path_i)
    else:
        # kill reaction path that not converge
        if os.path.isfile(f'{wd}/{path_i}/scan/xtbopt.xyz'):
            # # use the lowest point as end geom
            # lowgeom_l = lowestgeom(f'{wd}/{path_i}/scan/xtbopt.log')
            # with open(f'{wd}/{path_i}/coord.xyz', 'w') as f:
            #     f.seek(0, 0)
            #     for i in lowgeom_l:
            #         f.write(i)
            shutil.copy(f'{wd}/{path_i}/scan/xtbopt.xyz', f'{wd}/{path_i}/coord.xyz')
        elif os.path.isfile(f'{wd}/{path_i}/scan/NOT_CONVERGED'):
            print(f'{path_i}/scan not converging, check etemp setting')
            shutil.rmtree(f'{wd}/{path_i}', ignore_errors=True)
            del_l.append(path_i)
        else:
            print(f'{path_i}/scan failed')
            shutil.move(f'{wd}/{path_i}', f'{wd}/{path_i}_check')
            del_l.append(path_i)
    return del_l


def subprocess(wd, path_i, reactants, total_dict, rec, reactor, **kwargs):
    del_l = []
    xtb_setting = kwargs.get("xtb_setting")
    orca_setting = kwargs.get("orca_setting")
    homogene_dissoziation = kwargs.get("searching_setting").get('homogene', False)

    # seperate product for sucsses scan path
    if os.path.isdir(f'{wd}/{path_i}'):
        os.chdir(f'{wd}/{path_i}')
        seperate()

        fragment_l = [k for k in os.listdir(f'{wd}/{path_i}') if os.path.isdir(f'{wd}/{path_i}/{k}') and str.isdigit(k)]
        coord_old = pd.read_table(f'{wd}/{path_i}/scan/xtbopt.xyz', header=None, skiprows=2, sep='\s+')
        charg_file = pd.read_table(f'{wd}/{path_i}/scan/charges', header=None, sep='\s+')
        if not reactor:
            # update not reacting molecules
            reactants_tmp = total_dict[path_i]
            nproduct = len(fragment_l)
            count = 1
            for j in reactants:
                if j not in reactants_tmp and j in rec.keys():
                    rec[f'{path_i}/{str(count + nproduct)}'] = copy.deepcopy(rec[j])
                    rec[f'{path_i}/{str(count + nproduct)}']['name'] = \
                        rec[f'{path_i}/{str(count + nproduct)}']['mb'].split('|')[2]
                    # shutil.copytree(f'{wd}/{j}',f'{wd}/{path_i}/{str(count + nproduct)}')
                    count += 1

        # assign charg and m for reacted molecules
        uhf_d = get_uhf_from_excprep(f'{wd}/{path_i}/e.out', len(fragment_l))
        for j in fragment_l:
            rec[f'{path_i}/{j}'] = {}
            os.chdir(f'{wd}/{path_i}/{j}')
            coord_j = pd.read_table(f'{wd}/{path_i}/{j}/coord.xyz', header=None, skiprows=2, sep='\s+')
            nat_j = coord_j.shape[0]

            # assign charg
            charg_j = get_charge(coord_j, coord_old, charg_file)
            xyz_j = pd.read_table(f'{wd}/{path_i}/{j}/coord.xyz', header=None, skiprows=2, sep='\s+')

            # assign uhf
            uhf_j = uhf_d[int(j)]

            # xtb opt, try to calculate only with charge
            Struc(xyz_j, 1, charg_j).runxtb(f'{wd}/{path_i}/{j}/xtbopt', **xtb_setting)

            # correct uhf and charge (homogene_dissoziation)
            nel_j = get_nelectron(f'{wd}/{path_i}/{j}/xtbopt/xtb.out')
            if nel_j % 2 == 0 and uhf_j % 2 != 0:
                if homogene_dissoziation:
                    if charg_j > 0:
                        charg_j -= 1
                    else:
                        charg_j += 1
                else:
                    uhf_j -= 1
            elif nel_j % 2 != 0 and uhf_j % 2 == 0:
                if homogene_dissoziation:
                    if charg_j > 0:
                        charg_j -= 1
                    else:
                        charg_j += 1
                else:
                    uhf_j += 1

            # try to redo xtbopt
            rc = Struc(xyz_j, uhf_j + 1, charg_j).runxtb(f'{wd}/{path_i}/{j}/xtbopt', **xtb_setting)
            if rc != 0 :
                sys.exit(f'{wd}/{path_i}/{j} redo xtbopt wrong')
            rec = collect_sp_xyz(rec, f'{path_i}/{j}', f'{wd}/{path_i}/{j}/xtbopt', 'gfn2')

            if testcrop(wd, xyz_j):  # if struc stable
                # consider S and T for uhf = 0
                if uhf_j == 0:
                    if orca_setting.get('dft_level', 'r2scan-3c') == 'gfn2':
                        m_j = 0
                    else:
                        m_j = [1, 3]  # run both option to choose which one more stable
                else:
                    m_j = 2  # open shell structure, considerd as dublet

            else:
                m_j = 'inf'
                print(f'{wd}/{path_i}/{j}/xtbopt call xtb wrong')

            # update rec
            rec[f'{path_i}/{j}']['coord'] = xyz_j
            rec[f'{path_i}/{j}']['charge'] = int(charg_j)
            rec[f'{path_i}/{j}']['m'] = m_j  # hier m_j could be a list
            rec = collect_sp_xyz(rec, f'{path_i}/{j}', f'{wd}/{path_i}/{j}/xtbopt/', 'gfn2')

            # do molbar
            rec[f'{path_i}/{j}']['mb'], rec[f'{path_i}/{j}']['name'] = do_molbar(xyz_j, wd)

        # dft sp cal
        if orca_setting.get('dft_level', 'r2scan-3c') != 'gfn2':
            for j in fragment_l:
                os.chdir(f'{wd}/{path_i}/{j}')
                charg_j = rec[f'{path_i}/{j}']['charge']
                m_j = rec[f'{path_i}/{j}']['m']  # hier m could be a list
                xyz_j = rec[f'{path_i}/{j}']['coord']
                nat_j = xyz_j.shape[0]
                if isinstance(m_j, list):
                    for k in m_j:
                        if nat_j > 1 and k == 3:
                            # for m = 3 run xtbopt --uhf 2 again
                            Struc(xyz_j, 3, charg_j).runxtb(f'{wd}/{path_i}/{j}/orca/orca_{str(k)}', **xtb_setting)
                            rec = collect_sp_xyz(rec, f'{path_i}/{j}', f'{wd}/{path_i}/{j}/orca/orca_{str(k)}/', 'gfn2')

                        # run orca
                        Struc(xyz_j, k, charg_j).runorcasp(f'{wd}/{path_i}/{j}/orca/orca_{str(k)}',
                                                           caltyp='sp', **orca_setting)
                elif m_j != 'inf':
                    Struc(xyz_j, m_j, charg_j).runorcasp(f'{wd}/{path_i}/{j}/orca', caltyp='sp', **orca_setting)
                else:
                    del_l.append(path_i)
                    print(f"{path_i} failed due to m=inf (unstable product)")
            # collect dft results
            for j in fragment_l:
                tmp_l = ['orca/' + k for k in os.listdir(f'{wd}/{path_i}/{j}/orca')
                         if os.path.isdir(f'{wd}/{path_i}/{j}/orca/{k}') and "orca" in k]
                if not tmp_l:
                    tmp_l = ['orca']
                for k in tmp_l:
                    rec = collect_sp_xyz(rec, f'{path_i}/{j}', f'{wd}/{path_i}/{j}/{k}', 'dft')
    return del_l, rec



def elvauate_reaction(path_i, wd, rec, reactants_total_dict, reaction_rec, typ, **searching_setting):
    # reactants_total_dict: all reactaion paths in this cycle
    # del_l : reaction paths that not survior
    del_l = []
    reactor = searching_setting.get('reactor_like', False)
    energy_lim = searching_setting.get('elim', 300)
    if reactor:
        reactant_l = reactants_total_dict[path_i]
    else:
        reactant_l = [k for k in rec.keys() if k.split('/')[0] == reactants_total_dict[path_i][0].split('/')[0]]

    product_l = [k for k in rec.keys() if k.split('/')[0] == path_i]
    print(f"for path {path_i}: reactant = {reactant_l}; poduct = {product_l}")

    # generate reaction rec
    reactants_mb = [rec[k]['mb'] for k in reactant_l]
    product_mb = [rec[k]['mb'] for k in product_l]
    reaction_rec[path_i] = {}
    reaction_rec[path_i] = [reactants_mb, product_mb, 0]

    # test if reaction ends (products = reactants in one of cycles)
    if existreaction(reaction_rec, product_mb, 'r'):
        print(f"{path_i} failed due to reaction looped")
        del_l.append(f'{path_i}')
    elif existreaction(reaction_rec, product_mb, 'p'):
        print(f"{path_i} failed due to reaction exists")
        del_l.append(f'{path_i}')
    else:
        if len(product_l) - len(reactant_l) < 0 and \
                not all([testcrop(wd, rec[f'{k}']['coord'])
                         for k in product_l]):  # if unstable product (fall back to reactants)
            print(f"{path_i} failed due to unstable products")
            del_l.append(f'{path_i}')
        else:
            # filt with reaction energy
            ereactant = 0
            eproduct = 0
            # test with initial energy
            for i in [j for j in rec.keys() if len(j.split("/")[0]) == 1]:
                # for i in reactant_l:
                ereactant += rec[i]['sp_' + f'{typ}']
            for i in product_l:  # also filted out molecules carring with both charge and uhf
                eproduct += rec[f'{i}']['sp_' + f'{typ}']
                # if rec[f'{i}']['m'] != 1 and rec[f'{i}']['charge'] != 0:
                #     del_l.append(f'{path_i}')

            erec = eproduct - ereactant
            reaction_rec[path_i][2] = erec * 627.503
            print(f"{path_i}: reaction energy: {erec * 627.503}")
            if erec > energy_lim / 627.503:
                del_l.append(f'{path_i}')
                print(f"{path_i} failed due to recaction energy")
    return reaction_rec, del_l


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
