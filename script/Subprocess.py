import builtins
import os
import shutil
import sys
import time
import traceback
import copy
import pandas as pd
from read import existreaction
from get_molbar import do_molbar
from callsys import call_xtb, call_orca_dft
from operation import seperate, get_charge, creat_folder, write_xyz, \
    collect_energy, get_nelectron, testcrop, lowestgeom, get_uhf_from_excprep


def print(text):
    builtins.print(text)
    os.fsync(sys.stdout)


def subprocess(wd, path_i, reactants, total_dict, rec, method, alpb=None, reactor=False, mode=0):
    if alpb is not None:
        alpb = [alpb][0]

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
            print(f'{path_i}/scan not converging, check etmp setting')
            shutil.rmtree(f'{wd}/{path_i}', ignore_errors=True)
            del_l.append(path_i)
        else:
            print(f'{path_i}/scan failed')
            shutil.move(f'{wd}/{path_i}', f'{wd}/{path_i}_check')
            del_l.append(path_i)

    # seperate product for sucsses scan path
    if not del_l and os.path.isdir(f'{wd}/{path_i}'):
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
                    rec[f'{path_i}/{str(count + nproduct)}']['name'] = f'{path_i}/{str(count + nproduct)}'
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

            # xtb opt
            creat_folder(f'{wd}/{path_i}/{j}/xtbopt')  # not only for opt, also sp
            write_xyz(xyz_j, 'coord')
            call_xtb(f'{wd}/{path_i}/{j}/xtbopt', charg_j, nat_j, uhf=uhf_j, alpb=alpb, mode=mode)

            if os.path.isfile(f'{wd}/{path_i}/{j}/xtbopt/xtb.out'):
                nel = get_nelectron(f'{wd}/{path_i}/{j}/xtbopt/xtb.out')
                if nel % 2 == 0 and uhf_j % 2 != 0:  # even nel with open-shell
                    tmp_charg = charg_j
                    if len(fragment_l) ==1: # determined with charge equilibrium
                        charg_j = sum([rec[k]['charge'] for k in reactants])
                        uhf_j = 0
                    else: # excrop Priority
                        charg_j = 0
                    # try to redo xtbopt
                    call_xtb(f'{wd}/{path_i}/{j}/xtbopt', charg_j, nat_j, uhf=uhf_j, alpb=alpb, mode=mode)
                    if os.path.isfile(f'{wd}/{path_i}/{j}/xtbopt/xtb.out'):
                        nel = get_nelectron(f'{wd}/{path_i}/{j}/xtbopt/xtb.out')
                        if os.path.isfile(f'{wd}/{path_i}/{j}/xtbopt/xtbopt.xyz'):
                            xyz_j = pd.read_table(f'{wd}/{path_i}/{j}/xtbopt/xtbopt.xyz', header=None,
                                                  skiprows=2, sep='\s+')
                    else: # if reopt faild, roll back to old charge and uhf = 0
                        charg_j = tmp_charg
                        uhf_j = 0
                else:
                    if os.path.isfile(f'{wd}/{path_i}/{j}/xtbopt/xtbopt.xyz'):
                        xyz_j = pd.read_table(f'{wd}/{path_i}/{j}/xtbopt/xtbopt.xyz', header=None,
                                              skiprows=2, sep='\s+')
                if testcrop(wd, xyz_j): # if struc stable
                    m_j = uhf_j+1

                    #
                    #
                    # if nel % 2 != 0:
                    #     m_j = 2  # open shell structure, considerd as dublet
                    # else:
                    #     if method is not None:
                    #         m_j = [1, 3]  # run both option to choose which one more stable
                    #     else:
                    #         m_j = 1
                else:
                    m_j = 'inf'
            else:
                m_j = 'inf'
                print(f'{wd}/{path_i}/{j}/xtbopt call xtb wrong')

            # update rec
            rec[f'{path_i}/{j}']['coord'] = xyz_j
            rec[f'{path_i}/{j}']['charge'] = charg_j
            rec[f'{path_i}/{j}']['m'] = m_j  # hier m_j could be a list
            rec[f'{path_i}/{j}']['name'] = f'{path_i}/{j}'
            rec = collect_energy(rec, f'{wd}/{path_i}/{j}/xtbopt', f'{path_i}/{j}', "gfn2")

            # do molbar
            try:
                rec[f'{path_i}/{j}']['mb'] = do_molbar(xyz_j, wd)
            except:
                print("molbar error")
                traceback.print_exc()
                pass

        # dft cal
        for j in fragment_l:
            os.chdir(f'{wd}/{path_i}/{j}')
            charg = rec[f'{path_i}/{j}']['charge']
            m = rec[f'{path_i}/{j}']['m']  # hier m could be a list
            xyz_j = rec[f'{path_i}/{j}']['coord']
            nat_j = xyz_j.shape[0]
            if isinstance(m, list):
                for k in m:
                    creat_folder(f'{wd}/{path_i}/{j}/orca/orca_{str(k)}')
                    # for m = 3 run xtbopt --uhf 2 again
                    write_xyz(xyz_j, 'coord')
                    if nat_j > 1 and k == 3:
                        call_xtb(f'{wd}/{path_i}/{j}/orca/orca_{str(k)}', 1, nat_j, alpb=alpb, mode=mode,
                                 uhf=2)  # hier charged is no use
                        if os.path.isfile(f'{wd}/{path_i}/{j}/orca/orca_{str(k)}/xtbopt.xyz'):
                            xyz_j = pd.read_table(f'{wd}/{path_i}/{j}/xtbopt/xtbopt.xyz', header=None,
                                                  skiprows=2, sep='\s+')
                            write_xyz(xyz_j, 'coord')
                    print(f"doing dft for {path_i}/{j}/orca/orca_{str(k)}")
                    if method is not None:
                        call_orca_dft(f'{wd}/{path_i}/{j}/orca/orca_{str(k)}', method, charg, k, nat_j, alpb=alpb,
                                      mode=mode)
            elif m != 'inf':
                creat_folder(f'{wd}/{path_i}/{j}/orca')
                write_xyz(xyz_j, 'coord')
                print(f"doing dft for {path_i}/{j}")
                if method is not None:
                    call_orca_dft(f'{wd}/{path_i}/{j}/orca', method, charg, m, nat_j, alpb=alpb, mode=mode)
            else:
                del_l.append(path_i)
                print(f"{path_i} failed due to m=inf (unstable product)")
                del rec[f'{path_i}/{j}']
    return rec, del_l


def collect_dft(rec, wd):
    fragment_l = [k for k in list(rec.keys()) if "sp_dft" not in rec[k] and os.path.isdir(f'{wd}/{k}/orca')]
    if fragment_l:
        for j in fragment_l:
            os.chdir(f'{wd}/{j}')
            tmp_l = ['orca/' + k for k in os.listdir(f'{wd}/{j}/orca')
                     if os.path.isdir(f'{wd}/{j}/orca/{k}') and "orca" in k]
            if not tmp_l:
                tmp_l = ['orca']
            for k in tmp_l:
                if os.path.isfile(f'{wd}/{j}/{k}/orca.out'):
                    xyz_j = pd.read_table(f'{wd}/{j}/{k}/coord.xyz', header=None, skiprows=2, sep='\s+')
                    nat_j = xyz_j.shape[0]
                    if nat_j > 1:
                        if os.path.isfile(f'{wd}/{j}/{k}/orca.xyz'):
                            xyz_j = pd.read_table(f'{wd}/{j}/{k}/orca.xyz', header=None, skiprows=2, sep='\s+')
                    rec[f'{j}']['coord'] = xyz_j
                    # update rec and specify m and xyz
                    rec = collect_energy(rec, f'{wd}/{j}/{k}', f'{j}', "dft")
                else:
                    print("calling orca wrong")
    return rec


def elvauate_reaction(path_i, wd, rec, energy_lim, reactants_total_dict, reaction_rec, typ, reactor=False):
    # reactants_total_dict: all reactaion paths in this cycle
    # del_l : reaction paths that not survior
    del_l = []
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
