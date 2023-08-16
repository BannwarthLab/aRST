import copy
import os
import shutil
import time
import traceback
import copy
import billiard as multiprocessing
import pandas as pd
from read import existreaction
from get_molbar import do_molbar
# from callsys import seperate, get_charge, creat_folder, \
#      write_xyz, submit_xtb_opt, get_nelectron, submit_orca_sp, submit_orca_opt, submit_xtb_sp, collect_energy,testcrop

from submit_tmp import seperate, get_charge, creat_folder, \
   write_xyz, submit_xtb_opt, get_nelectron, submit_orca_sp, submit_orca_opt, submit_xtb_sp, collect_energy, testcrop


# def do_subprocess(wd, reactants, reactants_inter_dict, reactants_dis_dict, rec, method, energy_lim, alpb=None):
#     if alpb is not None:
#         alpb = [alpb][0]
#     # check scan stauts
#     reactants_inter_dict.update(reactants_dis_dict)  # now inter_dict = total dict
#     print(reactants, reactants_inter_dict)
#     path_l = list(reactants_inter_dict.keys())
#     print("before: ", path_l)
#     with multiprocessing.Pool(len(path_l)) as pool:
#         results = []
#         for i in path_l:
#             results.append(pool.apply_async(subprocess,
#                                             (wd, path_l, i, reactants, reactants_inter_dict, rec, method, alpb)))
#         pool.close()
#         pool.join()
#
#     del_l = []
#     for res in results:
#         rec.update(list(res.get())[0])
#         del_l += list(res.get())[1]
#
#     for i in list(set(del_l)):
#         path_l.remove(i)
#         del reactants_inter_dict[i]
#
#     print("after: ", path_l)
#     pool = multiprocessing.Pool(len(path_l))
#     results = [pool.apply_async(collect_dft, (rec, wd))]
#     pool.close()
#     pool.join()
#
#     for res in results:
#         rec.update(res.get())
#
#     try:
#         # get reaction energy
#         del_l = []
#         for path_i in path_l:
#             ereactant = 0
#             for i in reactants:
#                 ereactant += rec[i]['sp_dft']
#             eproduct = 0
#             product_l = [k for k in os.listdir(f'{wd}/{path_i}') if
#                          os.path.isdir(f'{wd}/{path_i}/{k}') and str.isdigit(k)]
#             for i in product_l:
#                 eproduct += rec[i]['sp_dft']
#
#             if eproduct - ereactant > energy_lim / 627.503:
#                 for i in product_l:
#                     print(f"rec energy too high del {path_i}/{i}")
#                     del_l.append(f'{path_i}/{i}')
#         for i in del_l:
#             if i in rec.keys():
#                 del rec[i]
#     except:
#         traceback.print_exc()
#         pass
#
#     return rec


def subprocess(wd, path_i, reactants, total_dict, rec, method, alpb=None, reactor=False):
    if alpb is not None:
        alpb = [alpb][0]
    print("i: ", path_i)
    del_l = []
    while not os.path.isfile(f'{wd}/{path_i}/scan/xtb.out'):
        time.sleep(20)
    else:
        # kill reaction path that not converge
        if os.path.isfile(f'{wd}/{path_i}/scan/xtbopt.xyz'):
            shutil.copy(f'{wd}/{path_i}/scan/xtbopt.xyz', f'{wd}/{path_i}/coord.xyz')
        elif os.path.isfile(f'{wd}/{path_i}/scan/NOT_CONVERGED'):
            shutil.rmtree(f'{wd}/{path_i}', ignore_errors=True)
            del_l.append(path_i)
        else:
            print(f'{path_i}/scan failed, maybe check electronic tmp')
            shutil.move(f'{wd}/{path_i}', f'{wd}/{path_i}_check')
            del_l.append(path_i)
    # seperate product for sucsses scan path
    if not del_l and os.path.isdir(f'{wd}/{path_i}'):
        os.chdir(f'{wd}/{path_i}')
        seperate()

        fragment_l = [k for k in os.listdir(f'{wd}/{path_i}') if os.path.isdir(f'{wd}/{path_i}/{k}') and str.isdigit(k)]
        print("fragment_l: ", fragment_l)
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
                    count += 1

        # assign charg and m for reacted molecules
        for j in fragment_l:
            rec[f'{path_i}/{j}'] = {}
            os.chdir(f'{wd}/{path_i}/{j}')
            coord_j = pd.read_table(f'{wd}/{path_i}/{j}/coord.xyz', header=None, skiprows=2, sep='\s+')
            nat_j = coord_j.shape[0]

            # assign charg
            charg_j = get_charge(coord_j, coord_old, charg_file)
            xyz_j = pd.read_table(f'{wd}/{path_i}/{j}/coord.xyz', header=None, skiprows=2, sep='\s+')

            # xtb opt
            creat_folder(f'{wd}/{path_i}/{j}/xtbopt')  # not only for opt, also sp
            write_xyz(xyz_j, 'coord')
            if nat_j > 1:
                submit_xtb_opt(f'{wd}/{path_i}/{j}/xtbopt', charg_j, alpb)
                while not os.path.isfile(f'{wd}/{path_i}/{j}/xtbopt/xtb.out'):
                    time.sleep(10)
                else:
                    rec = collect_energy(rec, f'{wd}/{path_i}/{j}/xtbopt', f'{path_i}/{j}', "gfn2")
                    if os.path.isfile(f'{wd}/{path_i}/{j}/xtbopt/xtbopt.xyz'):
                        xyz_j = pd.read_table(f'{wd}/{path_i}/{j}/xtbopt/xtbopt.xyz', header=None,
                                              skiprows=2, sep='\s+')

                # assgin m
                if get_nelectron(f'{wd}/{path_i}/{j}/xtbopt/xtb.out') % 2 != 0:  # open shell structure
                    m_j = 2
                else:
                    m_j = 1
            else:
                submit_xtb_sp(f'{wd}/{path_i}/{j}/xtbopt', charg_j, alpb)
                while not os.path.isfile(f'{wd}/{path_i}/{j}/xtbopt/xtb.out'):
                    time.sleep(2)
                else:
                    rec = collect_energy(rec, f'{wd}/{path_i}/{j}/xtbopt', f'{path_i}/{j}', "gfn2")

                if charg_j:
                    m_j = 1
                else:
                    m_j = 2

            # update rec
            rec[f'{path_i}/{j}']['coord'] = xyz_j
            rec[f'{path_i}/{j}']['charge'] = charg_j
            rec[f'{path_i}/{j}']['m'] = m_j
            rec[f'{path_i}/{j}']['name'] = f'{path_i}/{j}'

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
            m = rec[f'{path_i}/{j}']['m']
            xyz_j = rec[f'{path_i}/{j}']['coord']
            nat_j = xyz_j.shape[0]
            creat_folder(f'{wd}/{path_i}/{j}/orca')
            write_xyz(xyz_j, 'coord')
            if nat_j > 1:
                if charg != 0:
                    submit_orca_opt(f'{wd}/{path_i}/{j}/orca', method, charg, m, alpb)
                else:
                    submit_orca_opt(f'{wd}/{path_i}/{j}/orca', method, charg, m)
            else:
                if charg != 0:
                    submit_orca_sp(f'{wd}/{path_i}/{j}/orca', method, charg, m, alpb)
                else:
                    submit_orca_sp(f'{wd}/{path_i}/{j}/orca', method, charg, m)
    return rec, del_l


def collect_dft(rec, wd):
    fragment_l = [k for k in list(rec.keys()) if "sp_dft" not in rec[k] and os.path.isdir(f'{wd}/{k}/orca')]
    for j in fragment_l:
        os.chdir(f'{wd}/{j}')
        while not os.path.isfile(f'{wd}/{j}/orca/orca.out'):
            time.sleep(30)
            print(os.getcwd())
            print(f"i am hiere, {j}")
        else:
            xyz_j = pd.read_table(f'{wd}/{j}/orca/coord.xyz', header=None, skiprows=2, sep='\s+')
            nat_j = xyz_j.shape[0]
            if nat_j > 1:
                if os.path.isfile(f'{wd}/{j}/orca/orca.xyz'):
                    xyz_j = pd.read_table(f'{wd}/{j}/orca/orca.xyz', header=None, skiprows=2, sep='\s+')

            rec[f'{j}']['coord'] = xyz_j
            rec = collect_energy(rec, f'{wd}/{j}/orca', f'{j}', "dft")
    return rec


def elvauate_reaction(path_i, wd, rec, energy_lim, reactants_total_dict, reaction_rec,reactor=False):
    # reactants_total_dict: all reactaion paths in this cycle
    # del_l : reaction paths that not survior
    del_l = []
    if reactor:
        if type(reactants_total_dict[path_i]) == tuple:  # {'11': ('1/3', '1/5'), '12': ('1/3', '1/5'), '13': '1/1'}
            reactant_l = list(reactants_total_dict[path_i])
        else:
            reactant_l = [reactants_total_dict[path_i]]
    else:
        reactant_l = [k for k in rec.keys() if k.split('/')[0] == [reactants_total_dict[path_i]][0].split('/')[0]]

    product_l = [k for k in rec.keys() if k.split('/')[0] == path_i]

    # generate reaction rec
    reactants_mb = [rec[k]['mb'] for k in reactant_l]
    product_mb = [rec[k]['mb'] for k in product_l]
    reaction_rec[path_i] = {}
    reaction_rec[path_i] = [reactants_mb, product_mb, 0]

    # test if reaction ends (products = reactants in one of cycles)
    if existreaction(reaction_rec, product_mb):
        del_l.append(f'{path_i}')
    else:
        if len(product_l) - len(reactant_l) < 0 \
                and not testcrop(wd, rec[f'{path_i}/1']['coord']):  # unstable product (fall back to reactants)
            del_l.append(f'{path_i}')
        else:
            # filt with reaction energy
            ereactant = 0
            eproduct = 0
            for i in reactant_l:
                ereactant += rec[i]['sp_dft']
            for i in product_l:  # also filted out molecules carring with both charge and uhf
                eproduct += rec[f'{i}']['sp_dft']
                if rec[f'{i}']['m'] != 1 and rec[f'{i}']['charge'] != 0:
                    del_l.append(f'{path_i}')

            erec = eproduct - ereactant
            reaction_rec[path_i][2] = erec
            print(f"{path_i}: reaction energy: {erec * 627.503}")
            if erec > energy_lim / 627.503:
                del_l.append(f'{path_i}')

    return reaction_rec, del_l
