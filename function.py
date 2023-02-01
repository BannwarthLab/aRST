import sys
import time
import numpy as np
import pandas as pd
import os
import re


def get_xtbout_path(smi, wd_1, wd_2):
    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        path = f'{wd_2}/{smi}/xtb/ml/xtb.out'
    else:
        path = f'{wd_1}/{smi}/xtb/ml/xtb.out'
    return path


def get_coord_path(smi, wd_1, wd_2):
    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        path = f'{wd_2}/{smi}/xtb/ml/xtbopt.xyz'
    else:
        path = f'{wd_1}/{smi}/xtb/ml/xtbopt.xyz'
    return path


def get_ml_feature_path(smi, wd_1, wd_2):
    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        path = f'{wd_2}/{smi}/xtb/ml/ml_feature.csv'
    else:
        path = f'{wd_1}/{smi}/xtb/ml/ml_feature.csv'
    return path


def get_wbo_path(smi, wd_1, wd_2):
    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        path = f'{wd_2}/{smi}/xtb/ml/wbo'
    else:
        path = f'{wd_1}/{smi}/xtb/ml/wbo'
    return path


def get_uhf(smi, wd_2):
    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        smi_name = [item for item in filter(lambda x: x != '', smi.split('/'))][-1]
        with open(f'{wd_2}/{smi}/../e.out', 'r+') as f:
            l = []
            for line in f:
                if f"{smi_name}:" in line:
                    l.append(line)
                    # print(line)
            if len(l) >= 1:
                e = l[0]
            else:
                e = "0"
            if "UHF" in e:
                p1 = re.compile(r'[(](.*?)[)]', re.S)
                p2 = str(re.findall(p1, e))
                uhf = int(re.findall(r'\-?\d+\.?\d*', p2)[0])
            else:
                uhf = 0
        # if smi_name.find("disso") >= 0 or smi_name.find("push") >= 0:
        #     with open(f'{wd_2}/{smi}/sep/e.out', 'r+') as f:
        #         l = []
        #         for line in f:
        #             if f"{fragment}:" in line:
        #                 l.append(line)
        #                 print(line)
        #         if len(l) >= 1:
        #             e = l[0]
        #         else:
        #             e = "0"
        #         if "UHF" in e:
        #             p1 = re.compile(r'[(](.*?)[)]', re.S)
        #             p2 = str(re.findall(p1, e))
        #             uhf = int(re.findall(r'\-?\d+\.?\d*', p2)[0])
        #         else:
        #             uhf = 0
        #      # charge = float(pd.read_table(f'{wd_2}/{smi}/xtb/ml/charges', header=None).sum())
        # else:
        #     with open(f'{wd_2}/{smi}/../e.out', 'r+') as f:
        #         l = []
        #         for line in f:
        #             if f"{fragment}:" in line:
        #                 l.append(line)
        #                 print(line)
        #         if len(l)>=1:
        #             e = l[0]
        #         else:
        #             e = "0"
        #
        #         if "UHF" in e:
        #             p1 = re.compile(r'[(](.*?)[)]', re.S)
        #             p2 = str(re.findall(p1, e))
        #             uhf = int(re.findall(r'\-?\d+\.?\d*', p2)[0])
        #         else:
        #             uhf = 0
    else:
        uhf = 0
    # print(f"uhf={uhf}")
    # uhf = input("uhf = ")
    return uhf


def test_product(smi, wd_2):
    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        smi_name = [item for item in filter(lambda x: x != '', smi.split('/'))][-1]
        if smi_name.find("disso") >= 0 or smi_name.find("push") >= 0:
            with open(f'{wd_2}/{smi}/sep/e.out', 'r+') as f:
                l = []
                for line in f:
                    # line =''.join(line.split())
                    # if line.find("H2O") >= 0:
                    #     w = line.find("H2O")
                    #     if len(line) >= w+3:
                    #         if not line[w-1].isalpha() and not line[w+3].isalpha():
                    #             return 1
                    #         else:
                    #             return 0
                    #     else:
                    #         if not line[w-1].isalpha():
                    #             return 1
                    # else:
                    #     return 0

                    if line.find("C4H9CO2H") >= 0:
                        w = line.find("C4H9CO2H")
                        if len(line) >= w + 8:
                            if not line[w - 1].isalpha() and not line[w + 8].isalpha():
                                print("reaction completed")
                                sys.exit()
                            else:
                                return 0
                        else:
                            if not line[w - 8].isalpha():
                                print("reaction completed")
                                sys.exit()

        else:
            with open(f'{wd_2}/{smi}/../e.out', 'r+') as f:
                l = []
                for line in f:
                    line = ''.join(line.split())
                    # if line.find("H2O") >= 0:
                    #     w = line.find("H2O")
                    #     if len(line) >= w + 3:
                    #         if not line[w - 1].isalpha() and not line[w + 3].isalpha():
                    #             return 1
                    #         else:
                    #             return 0
                    #     else:
                    #         if not line[w - 1].isalpha():
                    #             return 1
                    # else:
                    #     return 0

                    if line.find("C4H9CO2H") >= 0:
                        w = line.find("C4H9CO2H")
                        if len(line) >= w + 8:
                            if not line[w - 1].isalpha() and not line[w + 8].isalpha():
                                print("reaction completed")
                                sys.exit()
                            else:
                                return 0
                        else:
                            if not line[w - 8].isalpha():
                                print("reaction completed")
                                sys.exit()
                    else:
                        return 0

    else:
        return 0


def get_xtbenergy(smi, wd_1, wd_2):
    with open(get_xtbout_path(smi, wd_1, wd_2), 'r') as f:
        l = []
        for line in f:
            if "TOTAL ENERGY" in line:
                l.append(line)
        e = l[0]
        singlepoint_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
    sp_kcal = singlepoint_eh * 627.503
    return sp_kcal


# def get_dftsp(smi, wd_1, wd_2):
#     if smi.find("disso") >= 0 or smi.find("push") >= 0:
#         with open(f'{wd_2}/{smi}/orca/r2scan-3c/orca.out', 'r') as f:
#             l = []
#             for line in f:
#                 if "FINAL SINGLE POINT ENERGY" in line:
#                     l.append(line)
#             e = l[0]
#             singlepoint_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
#         sp_kcal = singlepoint_eh * 627.503
#     else:
#         with open(f'{wd_1}/{smi}/orca/r2scan-3c/orca.out', 'r') as f:
#             l = []
#             for line in f:
#                 if "FINAL SINGLE POINT ENERGY" in line:
#                     l.append(line)
#             e = l[0]
#             singlepoint_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
#         sp_kcal = singlepoint_eh * 627.503
#     return sp_kcal


def get_dftsp(smi,wd_1,wd_2):

    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        with open(f'{wd_2}/{smi}/orca/r2scan-3c/orca.out', 'r') as f:
            l = []
            triger = False
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:
                    l.append(line)
            if not l:
                with open(f'{wd_2}/{smi}/orca/r2scan-3c/coord.xyz', 'r') as f:
                    for line in f:
                        if "H" in line:
                            singlepoint_eh = 0
            else:
                e = l[0]
                singlepoint_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
        sp_kcal = singlepoint_eh * 627.503
    else:
        with open(f'{wd_1}/{smi}/orca/r2scan-3c/orca.out', 'r') as f:
            l = []
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:
                    l.append(line)
            e = l[0]
            singlepoint_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
        sp_kcal = singlepoint_eh * 627.503
    return sp_kcal

def get_trv(smi,wd_1,wd_2):

    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        with open(f'{wd_2}/{smi}/xtb/gfn2xtb-fre/GTRV', 'r') as f:
            l = []
            for line in f:
                l.append(line)
            else:
                e = l[0]
                if 'E' in e or 'e' in e:
                    a = re.findall(r'(\-?\d+\.?\d*)([eE][-+]?[0-9]+)?$', e)[0]
                    trv_eh = float(a[0] + a[1])
                else:
                    trv_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])

        trv_kcal = trv_eh * 627.503
    else:
        with open(f'{wd_1}/{smi}/xtb/gfn2xtb-fre/GTRV', 'r') as f:
            l = []
            for line in f:
                l.append(line)
            e = l[0]
            if 'E' in e or 'e' in e:
                a = re.findall(r'(\-?\d+\.?\d*)([eE][-+]?[0-9]+)?$', e)[0]
                trv_eh = float(a[0] + a[1])
            else:
                trv_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
        trv_kcal = trv_eh * 627.503
    # print(smi, trv_eh)
    return trv_kcal


def get_solv(smi,wd_1,wd_2):

    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        with open(f'{wd_2}/{smi}/xtb/gfn2xtb-alpb-solv/gas.out', 'r') as f:
            l = []
            triger = False
            for line in f:
                if "TOTAL ENERGY" in line:
                    l.append(line)
            e = l[0]
            gas_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
        with open(f'{wd_2}/{smi}/xtb/gfn2xtb-alpb-solv/solv.out', 'r') as f:
            l = []
            triger = False
            for line in f:
                if "TOTAL ENERGY" in line:
                    l.append(line)
            e = l[0]
            sol_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
        sol_kcal = (gas_eh-sol_eh) * 627.503
    else:
        with open(f'{wd_1}/{smi}/xtb/gfn2xtb-alpb-solv/gas.out', 'r') as f:
            l = []
            for line in f:
                if "TOTAL ENERGY" in line:
                    l.append(line)
            e = l[0]
            gas_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
        with open(f'{wd_1}/{smi}/xtb/gfn2xtb-alpb-solv/solv.out', 'r') as f:
            l = []
            triger = False
            for line in f:
                if "TOTAL ENERGY" in line:
                    l.append(line)
            e = l[0]
            sol_eh = float(re.findall(r'\-?\d+\.?\d*', e)[0])
        sol_kcal = (gas_eh-sol_eh) * 627.503
    return sol_kcal


def get_freeenergy(smi,wd_1,wd_2):
    return get_dftsp(smi,wd_1,wd_2)+get_trv(smi,wd_1,wd_2)+get_solv(smi,wd_1,wd_2)

def spr_inter(reactants, wd_1, wd_2):
    # 1: hexyne	2: H2
    # {'2,1': [(9.61370222, [1, 15]), (9.73839884, [1, 14]), (15.28320733, [1, 13])]}
    if len(reactants) > 1:
        h = []
        l = []
        for smi in reactants:
            df = pd.read_csv(get_ml_feature_path(smi, wd_1, wd_2))
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


def spr_intra(smi, wd_1, wd_2):
    wbo_dict = {}
    with open(get_wbo_path(smi, wd_1, wd_2)) as f:
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
            df = pd.read_csv(f'{get_ml_feature_path(smi, wd_1, wd_2)}', index_col=0)
            homo = [np.array(df.loc[:, 'HOAO_a (eV)'].round(5))[atom_1 - 1],
                    np.array(df.loc[:, 'HOAO_b (eV)'].round(5))[atom_1 - 1],
                    np.array(df.loc[:, 'HOAO_a (eV)'].round(5))[atom_2 - 1],
                    np.array(df.loc[:, 'HOAO_b (eV)'].round(5))[atom_2 - 1]]
            lumo = [np.array(df.loc[:, 'LUAO_a (eV)'].round(5))[atom_1 - 1],
                    np.array(df.loc[:, 'LUAO_b (eV)'].round(5))[atom_1 - 1],
                    np.array(df.loc[:, 'LUAO_a (eV)'].round(5))[atom_2 - 1],
                    np.array(df.loc[:, 'LUAO_b (eV)'].round(5))[atom_2 - 1]]

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


def submit_scan(comand):
    # sbatch job to cluster, the submit.sh file is kept in autorec
    os.system(f'cp /home/guests/ychen/script/submit_script/submit_scan.sh .')
    with open('submit_scan.sh', 'r') as f:
        line = f.readlines()
    line[48] = f'{comand}'
    with open('submit_scan.sh', 'w') as f:
        for i in line:
            f.write(i)
    os.system('sbatch submit_scan.sh')


def submit_orca():
    os.system(f'cp /home/guests/ychen/script/submit_script/submit_orca.sh .')
    os.system('sbatch submit_orca.sh')


def r_distance(smi, nr, wd_1, wd_2):
    r0 = {'H': 0.32, 'D': 0.32, 'He': 0.46, 'Li': 1.20, 'Be': 0.94, 'B': 0.77, 'C': 0.75, 'N': 0.71, 'O': 0.63,
          'F': 0.64,
          'Ne': 0.67, 'Na': 1.40, 'Mg': 1.25, 'Al': 1.13, 'Si': 1.04, 'P': 1.10, 'S': 1.02, 'CL': 0.99, 'Ar': 0.96,
          'K': 1.76, 'Ca': 1.54, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11,
          'Ni': 1.10, 'Cu': 1.12, 'Zn': 1.18, 'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17,
          'Rb': 2.10, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
          'Pd': 1.20, 'Ag': 1.28, 'Cd': 1.36, 'In': 1.42, 'Sn': 1.40, 'Sb': 1.40, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31,
          'Cs': 2.32, 'Ba': 1.96, 'La': 1.80, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74, 'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68,
          'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.70, 'Lu': 1.62, 'Hf': 1.52,
          'Ta': 1.46, 'W': 0.95 * 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33,
          'Tl': 1.44,
          'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01, 'Ac': 1.86, 'Th': 1.75,
          'Pa': 1.69, 'U': 1.70, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65,
          'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57, 'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34,
          'Mt': 1.29, 'Ds': 1.21, 'Rg': 1.21, 'Cn': 1.22, 'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65,
          'Og': 1.57, 'Cl': 0.99}
    coord = pd.read_table(get_coord_path(smi, wd_1, wd_2), header=None, skiprows=2, sep='\s+')
    return r0[coord[0][int(nr) - 1]]


def get_principal_moment(coord, aktiv):
    r0 = {'H': 0.32, 'D': 0.32, 'He': 0.46, 'Li': 1.20, 'Be': 0.94, 'B': 0.77, 'C': 0.75, 'N': 0.71, 'O': 0.63,
          'F': 0.64,
          'Ne': 0.67, 'Na': 1.40, 'Mg': 1.25, 'Al': 1.13, 'Si': 1.04, 'P': 1.10, 'S': 1.02, 'CL': 0.99, 'Ar': 0.96,
          'K': 1.76, 'Ca': 1.54, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11,
          'Ni': 1.10, 'Cu': 1.12, 'Zn': 1.18, 'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17,
          'Rb': 2.10, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
          'Pd': 1.20, 'Ag': 1.28, 'Cd': 1.36, 'In': 1.42, 'Sn': 1.40, 'Sb': 1.40, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31,
          'Cs': 2.32, 'Ba': 1.96, 'La': 1.80, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74, 'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68,
          'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.70, 'Lu': 1.62, 'Hf': 1.52,
          'Ta': 1.46, 'W': 0.95 * 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33,
          'Tl': 1.44,
          'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01, 'Ac': 1.86, 'Th': 1.75,
          'Pa': 1.69, 'U': 1.70, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65,
          'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57, 'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34,
          'Mt': 1.29, 'Ds': 1.21, 'Rg': 1.21, 'Cn': 1.22, 'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65,
          'Og': 1.57, 'Cl': 0.99}
    t1, t2, t3, t4, t5, t6 = 0, 0, 0, 0, 0, 0
    # matrix for moments of inertia is of form, with CN in mass as weight
    # *           |   y**2+z**2                         |
    # *           |    -y*x       z**2+x**2             | -i =0
    # *           |    -z*x        -z*y       x**2+y**2 |
    cn = 0
    for i in range(coord.shape[0]):
        if i != aktiv:
            x = coord.iloc[i, 1] - coord.iloc[aktiv, 1]
            y = coord.iloc[i, 2] - coord.iloc[aktiv, 2]
            z = coord.iloc[i, 3] - coord.iloc[aktiv, 3]
            r = np.linalg.norm(np.array([x, y, z]))
            rcov_c = float(r0[f'{coord[0][aktiv]}'])
            rcov_i = float(r0[f'{coord[0][i]}'])
            rr = (rcov_i + rcov_c) / r
            cn = 1 / (1 + np.exp(- 16 * (rr - 1)))
            atmass = 0.5 * cn / (r ** 2)
            t1 = t1 + atmass * (y ** 2 + z ** 2)
            t2 = t2 + atmass * (-y * x)
            t3 = t3 + atmass * (z ** 2 + x ** 2)
            t4 = t4 + atmass * (-z * x)
            t5 = t5 + atmass * (-z * y)
            t6 = t6 + atmass * (x ** 2 + y ** 2)
    # print(cn)
    arep = np.zeros((3, 3))
    arep[0][0], arep[1][0], arep[1][1], arep[2][0], arep[2][1], arep[2][2] = t1, t2, t3, t4, t5, t6
    arep = arep + arep.T - np.diag(arep.diagonal())

    v, P = np.linalg.eigh(arep)
    # print(v)
    # print(P)
    return v, P


def angle(v1, v2):
    # clock
    norm = np.linalg.norm(v1) * np.linalg.norm(v2)
    if norm == 0:
        rho = 0
        t = 0
    else:
        rho = np.arcsin(np.cross(v1, v2) / norm)
        t = np.arccos(np.dot(v1, v2) / norm)
    if rho <= 0:
        return t
    else:
        return -t


def rot_vec(v, yaw, pitch):
    # rotate the vector
    df_z = np.array([[np.cos(yaw), np.sin(yaw), 0],
                     [-np.sin(yaw), np.cos(yaw), 0],
                     [0, 0, 1]], dtype=object)
    df_y = np.array([[np.cos(pitch), 0, -np.sin(pitch)],
                     [0, 1, 0],
                     [np.sin(pitch), 0, np.cos(pitch)]], dtype=object)
    oldcoord = np.array([v[0], v[1], v[2]])
    # M = np.dot(df_y,df_z)
    # newcoord= np.dot(M,oldcoord.T )
    newcoord = np.dot(df_y, np.dot(df_z, oldcoord.T))
    return newcoord


def rotate(df, yaw, pitch):
    # rotate the df
    # clockwise rotation yaw in the xy plane and clockwise rotation pitch in the xz plane
    for i in range(df.shape[0]):
        oldcoord = np.array([df.iloc[i, 1], df.iloc[i, 2], df.iloc[i, 3]])
        newcoord = rot_vec(oldcoord, yaw, pitch)
        df.iloc[i, 1] = newcoord[0]
        df.iloc[i, 2] = newcoord[1]
        df.iloc[i, 3] = newcoord[2]
    return df


def trans(df, aktiv, d):
    # translate the molecule with aktiv as center to (0,0,d)
    x = df.iloc[aktiv, 1]
    y = df.iloc[aktiv, 2]
    z = df.iloc[aktiv, 3] - d
    for i in range(df.shape[0]):
        df.iloc[i, 1] = df.iloc[i, 1] - x
        df.iloc[i, 2] = df.iloc[i, 2] - y
        df.iloc[i, 3] = df.iloc[i, 3] - z
    return df


def aligment(smi_a, smi_b, rec, d, wd_1, wd_2):
    coord_a = pd.read_table(get_coord_path(smi_a, wd_1, wd_2), header=None, skiprows=2, sep='\s+')
    coord_b = pd.read_table(get_coord_path(smi_b, wd_1, wd_2), header=None, skiprows=2, sep='\s+')

    coord_ta = trans(df=coord_a.copy(deep=True), aktiv=rec[0] - 1, d=0)
    coord_tb = trans(df=coord_b.copy(deep=True), aktiv=0, d=0)

    # v1, P1 = get_principal_moment(coord=coord_a, aktiv=rec[0] - 1)
    # v2, P2 = get_principal_moment(coord=coord_b, aktiv=rec[1] - 1)
    Iy_list_1l = []
    angle_list_1l = []
    Iy_list_1r = []
    angle_list_1r = []
    Iy_list_2l = []
    angle_list_2l = []
    Iy_list_2r = []
    angle_list_2r = []
    for yaw in np.linspace(0, 2 * np.pi, 30):
        for pitch in np.linspace(0, 2 * np.pi, 30):
            coord_ra = rotate(coord_ta.copy(), yaw, pitch)
            coord_rb = rotate(coord_tb.copy(), yaw, pitch)
            v1, P1 = get_principal_moment(coord=coord_ra, aktiv=rec[0] - 1)
            v2, P2 = get_principal_moment(coord=coord_rb, aktiv=rec[1] - 1)
            if np.round(np.cross(P1[0], P1[1])[2], 2) == 1:
                if np.round(P1[1][1], 2) == 1:
                    # left hand
                    Iy_list_1l.append(abs(np.cross(P1[0], P1[2])[1]))
                    angle_list_1l.append([yaw, pitch])
                else:
                    # right hand
                    Iy_list_1r.append(abs(np.cross(P1[0], P1[2])[1]))
                    angle_list_1r.append([yaw, pitch])
            if np.round(np.cross(P2[0], P2[1])[2], 2) == 1:
                if np.round(P2[1][1], 2) == 1:
                    # left hand
                    Iy_list_2l.append(abs(np.cross(P2[0], P2[2])[1]))
                    angle_list_2l.append([yaw, pitch])
                else:
                    # right hand
                    Iy_list_2r.append(abs(np.cross(P2[0], P2[2])[1]))
                    angle_list_2r.append([yaw, pitch])

    right_1 = False
    right_2 = False

    if Iy_list_1l:
        if not Iy_list_1r:
            yaw_1, pitch_1 = angle_list_1l[Iy_list_1l.index(max(Iy_list_1l))][0], \
                             angle_list_1l[Iy_list_1l.index(max(Iy_list_1l))][1]
        else:
            if min(Iy_list_1l) <= min(Iy_list_1r):
                right_1 = True
                yaw_1, pitch_1 = angle_list_1r[Iy_list_1r.index(max(Iy_list_1r))][0], \
                                 angle_list_1r[Iy_list_1r.index(max(Iy_list_1r))][1]
            else:
                yaw_1, pitch_1 = angle_list_1l[Iy_list_1l.index(max(Iy_list_1l))][0], \
                                 angle_list_1l[Iy_list_1l.index(max(Iy_list_1l))][1]
    else:
        if not Iy_list_1r:
            yaw_1, pitch_1 = 0, 0
        else:
            right_1 = True
            yaw_1, pitch_1 = angle_list_1r[Iy_list_1r.index(max(Iy_list_1r))][0], \
                             angle_list_1r[Iy_list_1r.index(max(Iy_list_1r))][1]

    if Iy_list_2l:
        if not Iy_list_2r:
            yaw_2, pitch_2 = angle_list_2l[Iy_list_2l.index(max(Iy_list_2l))][0], \
                             angle_list_2l[Iy_list_2l.index(max(Iy_list_2l))][1]
        else:
            if min(Iy_list_2l) <= min(Iy_list_2r):
                right_2 = True
                yaw_2, pitch_2 = angle_list_2r[Iy_list_2r.index(max(Iy_list_2r))][0], \
                                 angle_list_2r[Iy_list_2r.index(max(Iy_list_2r))][1]
            else:
                yaw_2, pitch_2 = angle_list_2l[Iy_list_2l.index(max(Iy_list_2l))][0], \
                                 angle_list_2l[Iy_list_2l.index(max(Iy_list_2l))][1]
    else:
        if not Iy_list_2r:
            yaw_2, pitch_2 = 0, 0
        else:
            right_2 = True
            yaw_2, pitch_2 = angle_list_2r[Iy_list_2r.index(max(Iy_list_2r))][0], \
                             angle_list_2r[Iy_list_2r.index(max(Iy_list_2r))][1]

    coord_ra = rotate(coord_ta.copy(), yaw_1, pitch_1)
    coord_rb = rotate(coord_tb.copy(), yaw_2, pitch_2)
    v1, P1 = get_principal_moment(coord=coord_ra, aktiv=rec[0] - 1)
    v2, P2 = get_principal_moment(coord=coord_rb, aktiv=rec[1] - 1)
    if right_1:
        coord_ra = rotate(coord_ra, angle([P1[1][0], P1[1][1]], [1, 0]), 0)
        P1 = rot_vec(P1, angle([P1[1][0], P1[1][1]], [1, 0]), 0)
    if right_2:
        coord_rb = rotate(coord_rb, angle([P2[1][0], P2[1][1]], [1, 0]), 0)
        P2 = rot_vec(P2, angle([P2[1][0], P2[1][1]], [1, 0]), 0)

    #
    coord_b_final = trans(coord_rb, aktiv=rec[1] - 1, d=d)
    product = pd.concat([coord_ra, coord_b_final])
    atomnr_1 = coord_a.shape[0]
    atmonr_2 = coord_b.shape[0]

    atomnr = str(atomnr_1 + atmonr_2)
    product.to_csv('product.xyz', sep='\t', index=False, header=None)
    with open('product.xyz', 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(atomnr + '\n' + '\n' + content)
    with open('cadm.bild', 'w+') as f2:
        f2.write('.color red\n')
        f2.write(f'.arrow 0 0 0 0 0 10 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 0 {P1[2][0]} {P1[2][1]} {P1[2][2]} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[2][0]} {P2[2][1]} {P2[2][2] + d} 0.05 0.150000\n')

        f2.write('.color green\n')
        f2.write(f'.arrow 0 0 0 10 0 0 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 0 {P1[0][0]} {P1[0][1]} {P1[0][2]} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[0][0]} {P2[0][1]} {P2[0][2] + d} 0.05 0.150000\n')

        f2.write('.color blue\n')
        f2.write(f'.arrow 0 0 0 0 10 0 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 0 {P1[1][0]} {P1[1][1]} {P1[1][2]} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[1][0]} {P2[1][1]} {P2[1][2] + d} 0.05 0.150000\n')
    return atomnr_1


def push_process(push_rec_p, distance_rec, uhf, c, gfnff):
    if len(push_rec_p) > 1:
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant=1\n'
                     f'$scan\n'
                     f'mode = concerted\n'
                     f'distance: {push_rec_p[0][0]},{push_rec_p[0][1]},{distance_rec[0]};{distance_rec[0]},{distance_rec[1]},30\n'
                     f'distance: {push_rec_p[1][0]},{push_rec_p[1][1]},{distance_rec[2]};{distance_rec[2]},{distance_rec[3]},30\n$end')
            # f'distance: {push_rec_p[1][0]},{push_rec_p[1][1]},{distance_rec[2]};{distance_rec[0]},{distance_rec[3]},30\n$end')
    else:
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant=1\n'
                     f'$scan\n'
                     f'distance: {push_rec_p[0][0]},{push_rec_p[0][1]},{distance_rec[0]};{distance_rec[0]},{distance_rec[1]},30\n$end')
    if gfnff == True:
        submit_scan(f'xtb product.xyz --opt --gfnff --uhf {uhf} -c {c}  --alpb h2o --input scan.inp > xtb.out\n')
    else:
        submit_scan(f'xtb product.xyz --opt --uhf {uhf} -c {c}  --alpb h2o --input scan.inp > xtb.out\n')
    time.sleep(1)
    return


def dissociate_process(smi_name, dis_rec, uhf, charge, wd_1, wd_2):
    os.system(f'cp {get_coord_path(smi_name, wd_1, wd_2)} ./coord.xyz')
    start_distance = atom_distance(smi_a=smi_name, smi_b=smi_name, nr_1=dis_rec[0], nr_2=dis_rec[1], wd_1=wd_1,
                                   wd_2=wd_2)
    with open('scan.inp', 'w') as f:
        f.write(f'$constrain\n'
                f'force constant=2\n'
                f'$scan\n'
                f'distance: {dis_rec[0]},{dis_rec[1]},{start_distance};{start_distance},{start_distance * 5},30\n$end')
    submit_scan(f'xtb coord.xyz --opt  --uhf {uhf} -c {charge} --alpb h2o --input scan.inp > xtb.out\n')
    time.sleep(1)


def excprep_process(smi, wd_2):
    if not os.path.isdir(f'{wd_2}/{smi}/sep'):
        os.makedirs(f'{wd_2}/{smi}/sep')
    os.chdir(f'{wd_2}/{smi}/sep')
    if os.path.isfile(f'{wd_2}/{smi}/scan/xtbopt.xyz'):
        os.system(f'cp {wd_2}/{smi}/scan/xtbopt.xyz ./coord.xyz')
    else:
        if os.path.isfile(f'{wd_2}/{smi}/scan/xtblast.xyz'):
            os.system(f'cp {wd_2}/{smi}/scan/xtblast.xyz ./coord.xyz')
        else:
            with open(f'{wd_2}/{smi}/scan/xtbscan.log', 'r') as f:
                line = f.readlines()
            atom = line[0]
            line.reverse()
            index_a = len(line) - line.index(atom) - 1
            line.reverse()
            with open(f'./coord.xyz', 'w+') as f:
                for i in line[index_a:len(line)]:
                    f.write(i)
    os.system(f'excprep coord.xyz > e.out')
    time.sleep(1)
    return


def ml_feature_process(smi, uhf, c, wd_1, wd_2):
    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        if not os.path.isdir(f'{wd_2}/{smi}/xtb/ml'):
            os.makedirs(f'{wd_2}/{smi}/xtb/ml')
        os.chdir(f'{wd_2}/{smi}/xtb/ml')
        os.system(f'cp {wd_2}/{smi}/coord.xyz ./coord.xyz')
        os.system(f'~/xtb_mloutput/_build/xtb coord.xyz --opt -c {c} --uhf {uhf}  --ml_feature>xtb.out')
    else:
        if not os.path.isdir(f'{wd_1}/{smi}/xtb/ml'):
            os.makedirs(f'{wd_1}/{smi}/xtb/ml')
        os.chdir(f'{wd_1}/{smi}/xtb/ml')
        os.system(f'cp {wd_1}/{smi}/struc.xyz ./coord.xyz')
        os.system(f'~/xtb_mloutput/_build/xtb coord.xyz --opt --uhf {uhf}  --ml_feature>xtb.out')
    time.sleep(1)
    return


def dft_opt(smi, c, m, wd_1, wd_2):
    if smi.find("disso") >= 0 or smi.find("push") >= 0:
        if not os.path.isdir(f'{wd_2}/{smi}/orca/r2scan-3c'):
            os.makedirs(f'{wd_2}/{smi}/orca/r2scan-3c')
        os.chdir(f'{wd_2}/{smi}/orca/r2scan-3c')
        with open('orca.inp', 'w') as f3:
            f3.write(f'! r2SCAN-3c sp OPT \n* xyzfile {c} {m} coord.xyz ')
        os.system(f'cp {wd_2}/{smi}/coord.xyz ./coord.xyz')
        submit_orca()
        while not os.path.isfile(f'{wd_2}/{smi}/orca/r2scan-3c/orca.out'):
            time.sleep(5)
    else:
        if not os.path.isdir(f'{wd_1}/{smi}/orca/r2scan-3c'):
            os.makedirs(f'{wd_1}/{smi}/orca/r2scan-3c')
            os.chdir(f'{wd_1}/{smi}/orca/r2scan-3c')
            with open('orca.inp', 'w') as f3:
                f3.write(f'! r2SCAN-3c sp OPT \n* xyzfile {c} {m} coord.xyz ')
            os.system(f'cp {wd_1}/{smi}/struc.xyz ./coord.xyz')
            submit_orca()
            while not os.path.isfile(f'{wd_1}/{smi}/orca/r2scan-3c/orca.out'):
                time.sleep(5)
        else:
            pass


def atom_distance(smi_a, smi_b, nr_1, nr_2, wd_1, wd_2):
    coord_1 = pd.read_table(get_coord_path(smi_a, wd_1, wd_2), header=None, skiprows=2, sep='\s+')
    coord_2 = pd.read_table(get_coord_path(smi_b, wd_1, wd_2), header=None, skiprows=2, sep='\s+')

    x = coord_1[1][int(nr_1) - 1] - coord_2[1][int(nr_2) - 1]
    y = coord_1[2][int(nr_1) - 1] - coord_2[2][int(nr_2) - 1]
    z = coord_1[3][int(nr_1) - 1] - coord_2[3][int(nr_2) - 1]
    distance = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    # if max([coord_1[1][int(nr_1)-1], coord_1[2][int(nr_1)-1], coord_1[3][int(nr_1)-1]]) < 1:
    #      distance = distance * 10

    return distance


def push(reactants_temp, push_rec, push_smi_a, push_smi_b, count, num, wd_1, wd_2):
    push_smi_a = reactants_temp[int(push_smi_a) - 1]
    push_smi_b = reactants_temp[int(push_smi_b) - 1]

    if not os.path.isdir(f'{wd_2}/{num}_push_{count}/scan'):
        os.makedirs(f'{wd_2}/{num}_push_{count}/scan')
    os.chdir(f'{wd_2}/{num}_push_{count}/scan')
    align_distance = np.round(
        (r_distance(push_smi_a, push_rec[0][0], wd_1, wd_2) + r_distance(push_smi_b, push_rec[0][1], wd_1, wd_2)) * 5,
        2)
    atomnnr_a = aligment(smi_a=push_smi_a, smi_b=push_smi_b, rec=push_rec[0], d=align_distance, wd_1=wd_1, wd_2=wd_2)
    if len(push_rec) > 1:
        push_rec_p = [[push_rec[0][0], push_rec[0][1] + atomnnr_a], [push_rec[1][0], push_rec[1][1] + atomnnr_a]]
        start_distance_1 = align_distance
        end_distance_1 = np.round(
            (r_distance(push_smi_a, push_rec[0][0], wd_1, wd_2) + r_distance(push_smi_b, push_rec[0][1], wd_1,
                                                                             wd_2)) * 0.9, 2)
        coord = pd.read_table(f'{wd_2}/{num}_push_{count}/scan/product.xyz', header=None, skiprows=2, sep='\s+')
        x = coord[1][push_rec[1][0] - 1] - coord[1][push_rec[1][1] - 1 + atomnnr_a]
        y = coord[2][push_rec[1][0] - 1] - coord[2][push_rec[1][1] - 1 + atomnnr_a]
        z = coord[3][push_rec[1][0] - 1] - coord[3][push_rec[1][1] - 1 + atomnnr_a]
        start_distance_2 = np.round(np.sqrt(x ** 2 + y ** 2 + z ** 2) * 1.1, 2)
        end_distance_2 = np.round(
            (r_distance(push_smi_a, push_rec[1][0], wd_1, wd_2) + r_distance(push_smi_b, push_rec[1][1], wd_1,
                                                                             wd_2)) * 0.9, 2)
        distance_rec = [start_distance_1, end_distance_1, start_distance_2, end_distance_2]
    else:
        push_rec_p = [[push_rec[0][0], push_rec[0][1] + atomnnr_a]]
        start_distance = align_distance
        end_distance = np.round(
            (r_distance(push_smi_a, push_rec[0][0], wd_1, wd_2) + r_distance(push_smi_b, push_rec[0][1], wd_1,
                                                                             wd_2)) * 0.85, 2)
        distance_rec = [start_distance, end_distance]
    c = 0
    uhf = 0
    for i in reactants_temp:
        uhf = uhf + get_uhf(i, wd_2)
    if uhf == 2:
        uhf = 0
    if uhf != 0:
        gfnff = True
    else:
        gfnff = False

    push_process(push_rec_p, distance_rec, uhf, c, gfnff)
    reactants_temp.remove(push_smi_a)
    reactants_temp.remove(push_smi_b)
    while not os.path.isfile(f'{wd_2}/{num}_push_{count}/scan/xtbscan.log'):
        time.sleep(5)

    excprep_process(f'{num}_push_{count}', wd_2)
    sp_kcal_temp = 0
    for i in os.listdir(f'{wd_2}/{num}_push_{count}/sep'):
        if os.path.isdir(f'{wd_2}/{num}_push_{count}/sep/{i}'):
            os.chdir(f'{wd_2}/{num}_push_{count}/sep/{i}')
            uhf = get_uhf(f'{num}_push_{count}/sep/{i}', wd_2)
            c = 0
            ml_feature_process(smi=f'{num}_push_{count}/sep/{i}', uhf=uhf, c=c, wd_1=wd_1, wd_2=wd_2)
            dft_opt(smi=f'{num}_push_{count}/sep/{i}', m=uhf + 1, wd_1=wd_1, wd_2=wd_2, c=0)
            if not os.path.isdir(f'{wd_2}/{num}_push_{count}/sep/{i}/xtb/gfn2xtb-alpb-solv'):
                os.makedirs(f'{wd_2}/{num}_push_{count}/sep/{i}/xtb/gfn2xtb-alpb-solv')
            os.chdir(f'{wd_2}/{num}_push_{count}/sep/{i}/xtb/gfn2xtb-alpb-solv')
            if os.path.isfile(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/orca.xyz'):
                os.system(f'cp {wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/orca.xyz ./coord.xyz')
            else:
                os.system(f'cp {wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/coord.xyz ./coord.xyz')

            os.system(f'xtb coord.xyz --gfn 2 --grad > gas.out')
            os.system(f'xtb coord.xyz --gfn 2 --alpb h2o  --grad > solv.out')
            time.sleep(1)
            if not os.path.isdir(f'{wd_2}/{num}_push_{count}/sep/{i}/xtb/gfn2xtb-fre'):
                os.makedirs(f'{wd_2}/{num}_push_{count}/sep/{i}/xtb/gfn2xtb-fre')
            os.chdir(f'{wd_2}/{num}_push_{count}/sep/{i}/xtb/gfn2xtb-fre')
            if os.path.isfile(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/orca.xyz'):
                os.system(f'cp {wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/orca.xyz ./coord.xyz')
            else:
                os.system(f'cp {wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/coord.xyz ./coord.xyz')
            os.system(f'xtb coord.xyz --gfn 2 --ohess > freq.out')
            os.system(f'TRVthermo --io xtb --rotthr 50 --fxyz coord.xyz --fvib freq.out --temp 298.15 > thermo.out')
            time.sleep(1)
            reactants_temp.append(f'{num}_push_{count}/sep/{i}')
            sp_kcal_temp = sp_kcal_temp + get_freeenergy(f'{num}_push_{count}/sep/{i}', wd_1, wd_2)
    # reactants_temp = list(i.replace(f'{num}', '') for i in reactants_temp)
    return reactants_temp, sp_kcal_temp


def dis(reactants_temp, dis_dict, count, num, wd_1, wd_2):
    # dis_gap = list(dis_dict[list(dis_dict.keys())[0]].keys())[0]
    dis_smi_name = reactants_temp[int(list(dis_dict.keys())[0]) - 1]
    dis_rec = dis_dict[list(dis_dict.keys())[0]][list(dis_dict[list(dis_dict.keys())[0]].keys())[0]]

    if not os.path.isdir(f'{wd_2}/{num}_disso_{count}/scan'):
        os.makedirs(f'{wd_2}/{num}_disso_{count}/scan')
    os.chdir(f'{wd_2}/{num}_disso_{count}/scan')
    uhf = get_uhf(f'{dis_smi_name}', wd_2)
    c = 0
    dissociate_process(dis_smi_name, dis_rec, uhf, c, wd_1, wd_2)
    time_1 = time.time()
    while not os.path.isfile(f'{wd_2}/{num}_disso_{count}/scan/xtbscan.log'):
        time.sleep(5)

    excprep_process(f'{num}_disso_{count}', wd_2)
    reactants_temp.remove(dis_smi_name)
    # reactants_temp = list(i.replace(f'{num}_', '') for i in reactants_temp)

    # Optimize fragments and do the ml_feature
    sp_kcal_temp = 0
    for i in os.listdir(f'{wd_2}/{num}_disso_{count}/sep'):
        if os.path.isdir(f'{wd_2}/{num}_disso_{count}/sep/{i}'):
            os.chdir(f'{wd_2}/{num}_disso_{count}/sep/{i}')
            uhf = get_uhf(f'{num}_disso_{count}/sep/{i}', wd_2)
            # uhf = 0
            c = 0
            dft_opt(smi=f'{num}_disso_{count}/sep/{i}', m=uhf + 1, wd_2=wd_2, wd_1=wd_1, c=0)
            while not os.path.isfile(f'{wd_2}/{num}_disso_{count}/sep/{i}/orca/r2scan-3c/orca.out'):
                time.sleep(5)
            ml_feature_process(smi=f'{num}_disso_{count}/sep/{i}', uhf=uhf, c=c, wd_1=wd_1, wd_2=wd_2)



            if not os.path.isdir(f'{wd_2}/{num}_disso_{count}/sep/{i}/xtb/gfn2xtb-alpb-solv'):
                os.makedirs(f'{wd_2}/{num}_disso_{count}/sep/{i}/xtb/gfn2xtb-alpb-solv')
            os.chdir(f'{wd_2}/{num}_disso_{count}/sep/{i}/xtb/gfn2xtb-alpb-solv')
            if os.path.isfile(f'{wd_2}/{num}_disso_{count}/sep/{i}/orca/r2scan-3c/orca.xyz'):
                os.system(f'cp {wd_2}/{num}_disso_{count}/sep/{i}/orca/r2scan-3c/orca.xyz ./coord.xyz')
            else:
                os.system(f'cp {wd_2}/{num}_disso_{count}/sep/{i}/orca/r2scan-3c/coord.xyz ./coord.xyz')
            os.system(f'xtb coord.xyz --gfn 2 --grad > gas.out')
            os.system(f'xtb coord.xyz --gfn 2 --alpb h2o  --grad > solv.out')
            time.sleep(1)
            if not os.path.isdir(f'{wd_2}/{num}_disso_{count}/sep/{i}/xtb/gfn2xtb-fre'):
                os.makedirs(f'{wd_2}/{num}_disso_{count}/sep/{i}/xtb/gfn2xtb-fre')
            os.chdir(f'{wd_2}/{num}_disso_{count}/sep/{i}/xtb/gfn2xtb-fre')
            if os.path.isfile(f'{wd_2}/{num}_disso_{count}/sep/{i}/orca/r2scan-3c/orca.xyz'):
                os.system(f'cp {wd_2}/{num}_disso_{count}/sep/{i}/orca/r2scan-3c/orca.xyz ./coord.xyz')
            else:
                os.system(f'cp {wd_2}/{num}_disso_{count}/sep/{i}/orca/r2scan-3c/coord.xyz ./coord.xyz')
            os.system(f'xtb coord.xyz --gfn 2 --ohess > freq.out')
            os.system(f'TRVthermo --io xtb --rotthr 50 --fxyz coord.xyz --fvib freq.out --temp 298.15 > thermo.out')
            time.sleep(1)


            reactants_temp.append(f'{num}_disso_{count}/sep/{i}')
            sp_kcal_temp = sp_kcal_temp + get_freeenergy(f'{num}_disso_{count}/sep/{i}', wd_1, wd_2)
    return reactants_temp, sp_kcal_temp


def push_BH76(reactants_temp, push_rec, push_smi_a, push_smi_b, count, num, wd_1, wd_2, c_l, m_l):
    push_smi_a = reactants_temp[int(push_smi_a) - 1]
    push_smi_b = reactants_temp[int(push_smi_b) - 1]

    if not os.path.isdir(f'{wd_2}/{num}_push_{count}/scan'):
        os.makedirs(f'{wd_2}/{num}_push_{count}/scan')
    os.chdir(f'{wd_2}/{num}_push_{count}/scan')
    align_distance = np.round(
        (r_distance(push_smi_a, push_rec[0][0], wd_1, wd_2) + r_distance(push_smi_b, push_rec[0][1], wd_1, wd_2)) * 5,
        2)
    atomnnr_a = aligment(smi_a=push_smi_a, smi_b=push_smi_b, rec=push_rec[0], d=align_distance, wd_1=wd_1, wd_2=wd_2)
    if len(push_rec) > 1:
        push_rec_p = [[push_rec[0][0], push_rec[0][1] + atomnnr_a], [push_rec[1][0], push_rec[1][1] + atomnnr_a]]
        start_distance_1 = align_distance
        end_distance_1 = np.round(
            (r_distance(push_smi_a, push_rec[0][0], wd_1, wd_2) + r_distance(push_smi_b, push_rec[0][1], wd_1,
                                                                             wd_2)) * 0.9, 2)
        coord = pd.read_table(f'{wd_2}/{num}_push_{count}/scan/product.xyz', header=None, skiprows=2, sep='\s+')
        x = coord[1][push_rec[1][0] - 1] - coord[1][push_rec[1][1] - 1 + atomnnr_a]
        y = coord[2][push_rec[1][0] - 1] - coord[2][push_rec[1][1] - 1 + atomnnr_a]
        z = coord[3][push_rec[1][0] - 1] - coord[3][push_rec[1][1] - 1 + atomnnr_a]
        start_distance_2 = np.round(np.sqrt(x ** 2 + y ** 2 + z ** 2) * 1.1, 2)
        end_distance_2 = np.round(
            (r_distance(push_smi_a, push_rec[1][0], wd_1, wd_2) + r_distance(push_smi_b, push_rec[1][1], wd_1,
                                                                             wd_2)) * 0.9, 2)
        distance_rec = [start_distance_1, end_distance_1, start_distance_2, end_distance_2]
    else:
        push_rec_p = [[push_rec[0][0], push_rec[0][1] + atomnnr_a]]
        start_distance = align_distance
        end_distance = np.round(
            (r_distance(push_smi_a, push_rec[0][0], wd_1, wd_2) + r_distance(push_smi_b, push_rec[0][1], wd_1,
                                                                             wd_2)) * 0.9, 2)
        distance_rec = [start_distance, end_distance]
    c = 0
    for k in c_l:
        c = int(c + k)
    uhf = 0
    for i in reactants_temp:
        uhf = uhf + get_uhf(i, wd_2)
    if uhf == 2:
        uhf = 0
    if uhf != 0:
        gfnff = True
    else:
        gfnff = False

    push_process(push_rec_p, distance_rec, uhf, c, gfnff)
    reactants_temp.remove(push_smi_a)
    reactants_temp.remove(push_smi_b)
    while not os.path.isfile(f'{wd_2}/{num}_push_{count}/scan/xtbscan.log'):
        time.sleep(5)

    excprep_process(f'{num}_push_{count}', wd_2)
    sp_kcal_temp = 0
    for i in os.listdir(f'{wd_2}/{num}_push_{count}/sep'):
        if os.path.isdir(f'{wd_2}/{num}_push_{count}/sep/{i}'):
            os.chdir(f'{wd_2}/{num}_push_{count}/sep/{i}')
            uhf = get_uhf(f'{num}_push_{count}/sep/{i}', wd_2)
            # c = 0
            # ml_feature_process(smi=f'{num}_push_{count}/sep/{i}', uhf=uhf, c=c, wd_1=wd_1, wd_2=wd_2)
            # dft_opt(smi=f'{num}_push_{count}/sep/{i}', m=uhf +1, wd_1=wd_1, wd_2=wd_2,c = 0)
            # reactants_temp.append(f'{num}_push_{count}/sep/{i}')
            # sp_kcal_temp = sp_kcal_temp + get_dftsp(f'{num}_push_{count}/sep/{i}', wd_1, wd_2)

            m = m_l[0]
            c = c_l[0]
            if not os.path.isdir(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c'):
                os.makedirs(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c')
            os.chdir(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c')
            with open('orca.inp', 'w') as f3:
                f3.write(f'! r2SCAN-3c sp  \n* xyzfile {c} {m} coord.xyz ')
            os.system(f'cp {wd_2}/{num}_push_{count}/sep/{i}/coord.xyz ./coord.xyz')
            submit_orca()
            time_1 = time.time()
            while not os.path.isfile(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/orca.out'):
                time.sleep(5)
                time_2 = time.time()
                if time_2 - time_1 > 100:
                    pass

            with open(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/orca.out', 'r') as f:
                l = []
                for line in f:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        l.append(line)
            if not l:
                m = m_l[1]
                c = c_l[1]
                os.chdir(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c')
                os.system('rm orca.out')
                os.system('rm orca.inp')
                with open('orca.inp', 'w') as f3:
                    f3.write(f'! r2SCAN-3c sp  \n* xyzfile {c} {m} coord.xyz ')
                os.system(f'cp {wd_2}/{num}_push_{count}/sep/{i}/coord.xyz ./coord.xyz')
                submit_orca()
                time_1 = time.time()
                while not os.path.isfile(f'{wd_2}/{num}_push_{count}/sep/{i}/orca/r2scan-3c/orca.out'):
                    time.sleep(15)
                    time_2 = time.time()
                    if time_2 - time_1 > 100:
                        pass

            ml_feature_process(smi=f'{num}_push_{count}/sep/{i}', uhf=uhf, c=c, wd_1=wd_1, wd_2=wd_2)
            reactants_temp.append(f'{num}_push_{count}/sep/{i}')
            sp_kcal_temp = sp_kcal_temp + get_dftsp(f'{num}_push_{count}/sep/{i}', wd_1, wd_2)

    # reactants_temp = list(i.replace(f'{num}', '') for i in reactants_temp)
    return reactants_temp, sp_kcal_temp
