import os
import shutil

import numpy as np
import pandas as pd

"""Please place this script and your .xyz file in the same folder

please load module : xtb/ml-hessian !!!!!


if you want to predict the intermolecular reaction, fill in the target_name_1/2 
e.g. target_name_1 = str("ch3br"), target_name_2 = str("cl");
if you want to predict the intramolecular reaction, just fill in the target_name_1 
e.g. target_name_1 = str("ch3br"), target_name_2 = str("").
"""

target_name_1 = str("a")
c_1 = 0
u_1 = 0  # !!! please change here if charge and uhf are not 0!
target_name_2 = str("b")
c_2 = 0
u_2 = 0  # !!! please change here if charge and uhf are not 0!


#######################Functions#########################
def spr_inter(reactants, wd):
    # 1: hexyne	2: H2
    # {'2,1': [(9.61370222, [1, 15]), (9.73839884, [1, 14]), (15.28320733, [1, 13])]}
    if len(reactants) > 1:
        h = []
        l = []
        for smi in reactants:
            df = pd.read_csv(f'{wd}/{smi}_ml/ml_feature.csv')
            if 'HOAO_a (eV)' in df.columns:
                h.append(df['HOAO_a (eV)'].tolist())
                h.append(df['HOAO_b (eV)'].tolist())
                l.append(df['LUAO_a (eV)'].tolist())
                l.append(df['LUAO_b (eV)'].tolist())
            else:
                h.append(df['HOAO (eV)'].tolist())
                h.append(df['HOAO (eV)'].tolist())
                l.append(df['LUAO (eV)'].tolist())
                l.append(df['LUAO (eV)'].tolist())
        d = {}
        d_1 = {}
        d_2 = {}
        for i in range(len(h)):
            if i % 2 == 0:
                start = 1
            else:
                start = 0
            for j in range(*{'start': start, 'stop': len(l), 'step': 2}.values()):
                if int(i / 2) != int(j / 2):
                    h_i = h[i]
                    l_j = l[j]
                    # print(h_i)
                    # print(l_j)
                    d[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = {}
                    d_2[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = {}
                    for one in range(len(h_i)):
                        for two in range(len(l_j)):
                            try:
                                h_i[one] = float(h_i[one])
                            except ValueError:
                                h_i[one] = float('inf')
                            try:
                                l_j[two] = float(l_j[two])
                            except ValueError:
                                l_j[two] = float('inf')

                            key = np.round(abs(float(h_i[one]) - float(l_j[two])), 3)
                            if key in d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].keys():
                                d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key].append([one + 1, two + 1])
                            else:
                                d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key] = [[one + 1, two + 1]]

                            for item_2 in list(d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].keys()):
                                if 0 < abs(key - item_2) < 0.01:
                                    temp = list(d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][item_2])
                                    temp[0].append('*')
                                    d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key] = \
                                        d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][
                                            key] + temp
                                    del d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][item_2]
                    # print(d[f'{int(i/2)+1},{int(j/2)+1}'])
                    # print(int(i/2)+1, int(j/2)+1)
                    # Discard excess data and organize it
                    if f'{int(i / 2) + 1},{int(j / 2) + 1}' in d_1.keys():
                        d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] + sorted(
                            d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].items())
                    else:
                        d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = sorted(
                            d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].items())
                    # print(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'])
                    if len(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}']) > 1:
                        if d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][0][0] < \
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][1][
                                    0]:
                            d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'].remove(
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][1])
                        else:
                            d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'].remove(
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][0])
                    for u in list(dict(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}']).keys()):
                        l_1 = []
                        l_2 = []
                        for item in list(dict(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'])[u]):
                            if (item[0] not in l_2) and (item[1] not in l_2):
                                l_2.append(item[0])
                                l_2.append(item[1])
                                l_1.append(item)
                                d_2[f'{int(i / 2) + 1},{int(j / 2) + 1}'][u] = l_1
                                # dict(d_1[f'{i + 1},{j + 1}'])[u] = l_1
        for a in d_2.keys():
            for b in list(d_2[a].keys())[2:]:
                del d_2[a][b]
        # print("push_dict: alpha beta: ")
        # print(d_2)

        d = {}
        d_1 = {}
        d_3 = {}
        for i in range(len(h)):
            if i % 2 != 0:
                start = 1
            else:
                start = 0
            for j in range(*{'start': start, 'stop': len(l), 'step': 2}.values()):
                if int(i / 2) != int(j / 2):
                    h_i = h[i]
                    l_j = l[j]
                    # print(h_i)
                    # print(l_j)
                    d[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = {}
                    d_3[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = {}
                    for one in range(len(h_i)):
                        for two in range(len(l_j)):
                            try:
                                h_i[one] = float(h_i[one])
                            except ValueError:
                                h_i[one] = float('inf')
                            try:
                                l_j[two] = float(l_j[two])
                            except ValueError:
                                l_j[two] = float('inf')
                            key = np.round(abs(float(h_i[one]) - float(l_j[two])), 3)
                            if key in d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].keys():
                                d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key].append([one + 1, two + 1])
                            else:
                                d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key] = [[one + 1, two + 1]]

                            for item_2 in list(d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].keys()):
                                if 0 < abs(key - item_2) < 0.01:
                                    temp = list(d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][item_2])
                                    temp[0].append('*')
                                    d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key] = \
                                        d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][
                                            key] + temp
                                    del d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][item_2]

                    if f'{int(i / 2) + 1},{int(j / 2) + 1}' in d_1.keys():
                        d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] + sorted(
                            d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].items())
                    else:
                        d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = sorted(
                            d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].items())
                    # print(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'])
                    if len(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}']) > 1:
                        if d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][0][0] < \
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][1][
                                    0]:
                            d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'].remove(
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][1])
                        else:
                            d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'].remove(
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][0])
                    for u in list(dict(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}']).keys()):
                        l_1 = []
                        l_2 = []
                        for item in list(dict(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'])[u]):
                            if (item[0] not in l_2) and (item[1] not in l_2):
                                l_2.append(item[0])
                                l_2.append(item[1])
                                l_1.append(item)
                                d_3[f'{int(i / 2) + 1},{int(j / 2) + 1}'][u] = l_1
                                # dict(d_1[f'{i + 1},{j + 1}'])[u] = l_1
        for a in d_3.keys():
            for b in list(d_3[a].keys())[2:]:
                del d_3[a][b]
        # print("push_dict: alpha alpha: ")
        # print(d_3)

        return d_2, d_3
    else:
        d_2 = {}
        d_3 = {}
        return d_2, d_3


def spr_intra(smi, wd):
    wbo_dict = {}
    with open(f'{wd}/{smi}_ml/wbo') as f:
        lines = f.readlines()
        if len(lines) > 0:
            for item in lines:
                a, b, c = str.split(item)
                wbo_dict[c] = a, b
    order = list(wbo_dict.keys())
    gap = float('inf')
    atom_1_final, atom_2_final = float('inf'), float('inf')
    for item in order:
        if float(item) < 1.1:
            atom_1, atom_2 = int(list(wbo_dict[item])[0]), int(list(wbo_dict[item])[1])
            df = pd.read_csv(f'{wd}/{smi}_ml/ml_feature.csv', index_col=0)
            if 'HOAO_a (eV)' in df.columns:
                homo = [np.array(df.loc[:, 'HOAO_a (eV)'].round(5))[atom_1 - 1],
                        np.array(df.loc[:, 'HOAO_b (eV)'].round(5))[atom_1 - 1],
                        np.array(df.loc[:, 'HOAO_a (eV)'].round(5))[atom_2 - 1],
                        np.array(df.loc[:, 'HOAO_b (eV)'].round(5))[atom_2 - 1]]
                lumo = [np.array(df.loc[:, 'LUAO_a (eV)'].round(5))[atom_1 - 1],
                        np.array(df.loc[:, 'LUAO_b (eV)'].round(5))[atom_1 - 1],
                        np.array(df.loc[:, 'LUAO_a (eV)'].round(5))[atom_2 - 1],
                        np.array(df.loc[:, 'LUAO_b (eV)'].round(5))[atom_2 - 1]]
            else:
                homo = [np.array(df.loc[:, 'HOAO (eV)'].round(5))[atom_1 - 1],
                        np.array(df.loc[:, 'HOAO (eV)'].round(5))[atom_1 - 1],
                        np.array(df.loc[:, 'HOAO (eV)'].round(5))[atom_2 - 1],
                        np.array(df.loc[:, 'HOAO (eV)'].round(5))[atom_2 - 1]]
                lumo = [np.array(df.loc[:, 'LUAO (eV)'].round(5))[atom_1 - 1],
                        np.array(df.loc[:, 'LUAO (eV)'].round(5))[atom_1 - 1],
                        np.array(df.loc[:, 'LUAO (eV)'].round(5))[atom_2 - 1],
                        np.array(df.loc[:, 'LUAO (eV)'].round(5))[atom_2 - 1]]
            list_temp_1 = [abs(homo[0] - lumo[3]), abs(homo[1] - lumo[2])]
            gap_1 = min(list_temp_1)
            list_temp_2 = [abs(homo[2] - lumo[1]), abs(homo[3] - lumo[0])]
            gap_2 = min(list_temp_2)
            if gap_1 <= gap_2:
                if gap_1 <= gap:
                    atom_1_final, atom_2_final = atom_1, atom_2
                    gap = gap_1
            else:
                if gap_2 <= gap:
                    atom_1_final, atom_2_final = atom_1, atom_2
                    gap = gap_2
    # print(atom_1_final, atom_2_final, gap)

    return gap, [atom_1_final, atom_2_final]


################################################

wd = os.getcwd()
d = {}
d[target_name_1] = [c_1, u_1]
if target_name_2:
    d[target_name_2] = [c_2, u_2]

# call xTB
# os.system('module load xtb/ml-hessian')

for i in list(d.keys()):
    if os.path.isdir(f'{wd}/{i}_ml'):
        shutil.rmtree(f'{wd}/{i}_ml')
    os.makedirs(f'{wd}/{i}_ml')
    os.system(f'cp {wd}/{i}.xyz {wd}/{i}_ml')
    os.chdir(f'{wd}/{i}_ml')
    os.system(f'xtb {i}.xyz -c {d[i][0]} --uhf {d[i][1]} --ml_feature>xtb.out')

# output
if target_name_2:
    print('inter: ')
    print(spr_inter(list(d.keys()), wd))
    # e.g.: {'1,2': {11.003: [[2, 1]], 18.058: [[2, 2]]}, '2,1': {9.774000000000001: [[2, 2]], 11.34: [[1, 1]]}}
    # '1,2': first molecule provides HOAO and second molecule provides LUAO
    # the lowest HL_gap = 11.003 eV,
    # between 2. atom from first molecule and 1. atom from second molecule (atom index of original .xyz file),
    # since both molecule is exprected to provide HOAO/LUAO, so please check both options, that is [2,1] and [2,2].
    for k in list(d.keys()):
        print(f'intra_{k}: ')
        print(spr_intra(k, wd))
        # (9.16089, [1, 2]) : the lowest HL_gap = 9.16089 eV, between 1. and 2. atom

