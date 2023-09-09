import os
import traceback
import copy
import numpy as np
import pandas as pd
import ast
import random
import math


def find_value(dict, fvalue):
    result = False
    for key in dict.keys():
        if dict[key] == fvalue:
            result = True
            break
    return result


def get_key(dict, value):
    return [k for k, v in dict.items() if v == value][0]


def litgap(react_d, hoao_df, luao_df):
    tmp1 = luao_df.min().tolist()
    tmp2 = hoao_df.max().tolist()
    tmp2c = copy.deepcopy(tmp2)
    tmp2c.pop(tmp1.index(min(tmp1)))
    lit1 = min(tmp1) - max(tmp2c)

    tmp1c = copy.deepcopy(tmp1)
    tmp1c.pop(tmp2.index(max(tmp2)))
    lit2 = min(tmp1c) - max(tmp2c)
    del_l = [i for i in react_d.keys() if i < min(lit1, lit2)]
    for i in del_l:
        del react_d[i]
    return react_d
def build_gap_d(react_d, h_df, l_df):
    # {6.826: ['[1, 2]', '[[1, 1]]'], 8.268: ['[1, 2]', '[[3, 1]]'], 7.477: ['[2, 1]', '[[1, 2]]'], 7.878: ['[1, 2]',
    # '[[1, 2]]']}

    for i in list(h_df.columns):
        for j in list(l_df.columns):
            if i != j:
                struc1 = i + 1
                struc2 = j + 1
                for h_index in h_df.iterrows():
                    for l_index in l_df.iterrows():
                        at1 = h_index[0]
                        at2 = l_index[0]
                        a = [1]

                        e1 = float(h_df.loc[at1][i])
                        e2 = float(l_df.loc[at2][j])
                        gap = np.round(e2 - e1,4)
                        at1 += 1  # index to atom number
                        at2 += 1
                        value1 = [f'{[struc1, struc2]}', f'{[[at1, at2]]}']
                        value2 = [f'{[struc2, struc1]}', f'{[[at2, at1]]}']
                        try:
                            if gap not in react_d.keys():
                                # delete dupilicated pairs
                                if not find_value(react_d, value1) and not find_value(
                                        react_d, value2) and not np.isnan(gap):
                                    react_d[gap] = value1
                                elif find_value(react_d, value1) and get_key(react_d, value1) > gap:
                                    del react_d[get_key(react_d, value1)]
                                    react_d[gap] = value1
                                elif find_value(react_d, value2) and get_key(react_d, value2) > gap:
                                    del react_d[get_key(react_d, value2)]
                                    react_d[gap] = value2
                            else:  # equivalent pairs
                                if not np.isnan(gap):
                                    gap += random.uniform(1e-11, 1e-10)
                                    react_d[gap] = value1

                        except AttributeError:
                            react_d = {gap: value1}

    # print(react_d)
    # merge equivalent pairs
    result_dict = {}
    for gap in react_d.keys():
        value = react_d[gap]
        # print(f'gap = {gap}, value = {value}')
        gap_appoxi = np.round(gap, 4)
        if gap_appoxi not in result_dict.keys():
            result_dict[gap_appoxi] = value
        else:
            value_old = result_dict[gap_appoxi]
            at_l = ast.literal_eval(value_old[1])
            struc1_old = ast.literal_eval(value_old[0])[0]
            struc2_old = ast.literal_eval(value_old[0])[1]
            struc1 = ast.literal_eval(value[0])[0]
            struc2 = ast.literal_eval(value[0])[1]
            at1 = ast.literal_eval(value[1])[0][0]
            at2 = ast.literal_eval(value[1])[0][1]
            for item in at_l:
                if at1 not in item and at2 not in item and [at1, at2] not in at_l and [at2, at1] not in at_l:
                    if struc1 == struc1_old and struc2 == struc2_old:
                        at_l.append([at1, at2])
                        value_old[1] = str(at_l)
                        result_dict[gap_appoxi] = value_old
                    elif struc2 == struc1_old and struc1 == struc2_old:
                        at_l.append([at2, at1])
                        value_old[1] = str(at_l)
                        result_dict[gap_appoxi] = value_old
    return result_dict


def spr_inter_new(reactants, rec):
    hoaoa_d = {}
    at_l1 = []
    hoaob_d = {}
    at_l2 = []

    luaoa_d = {}
    at_l3 = []
    luaob_d = {}
    at_l4 = []
    for i in range(len(reactants)):
        hoaoa_d[i] = {}
        hoaob_d[i] = {}
        name = reactants[i]
        # if '_' in smi:
        #     smi = smi.split('_')[0]
        ml_feature = rec[name]['ml_feature']
        if 'hoao_a' in list(ml_feature.columns):
            hoaoa_l = list(ml_feature['hoao_a'])
            hoaob_l = list(ml_feature['hoao_b'])

        else:
            hoaoa_l = list(ml_feature['hoao'])
            hoaob_l = list(ml_feature['hoao'])
        max1 = max(hoaoa_l)
        target_hoao_a = list(reversed(sorted([k for k in hoaoa_l if max1 - k < 2.5])))
        max2 = max(hoaob_l)
        target_hoao_b = list(reversed(sorted([k for k in hoaob_l if max2 - k < 2.5])))
        for x, y in zip(target_hoao_a, target_hoao_b):
            if hoaoa_l.count(x) > 1:
                at_x = -1
                for c in range(0, hoaoa_l.count(x)):
                    at_x = hoaoa_l.index(x, at_x + 1)
                    hoaoa_d[i][at_x] = x
                    at_l1.append(at_x)
            else:
                at_x = hoaoa_l.index(x)
                hoaoa_d[i][at_x] = x
                at_l1.append(at_x)
            if hoaob_l.count(y) > 1:
                at_y = -1
                for c in range(0, hoaob_l.count(y)):
                    at_y = hoaob_l.index(y, at_y + 1)
                    hoaob_d[i][at_y] = x
                    at_l2.append(at_y)
            else:
                at_y = hoaob_l.index(y)
                hoaob_d[i][at_y] = y
                at_l2.append(at_y)
            at_l1 = list(set(at_l1))
            at_l2 = list(set(at_l2))

        luaoa_d[i] = {}
        luaob_d[i] = {}
        if 'luao_a' in ml_feature.columns:
            luaoa_l = list(ml_feature['luao_a'])
            luaob_l = list(ml_feature['luao_b'])

        else:
            luaoa_l = list(ml_feature['luao'])
            luaob_l = list(ml_feature['luao'])
        min1 = min(luaoa_l)
        target_luao_a = list(sorted([k for k in luaoa_l if k - min1 < 2.5]))
        min2 = min(luaob_l)
        target_luao_b = list(sorted([k for k in luaob_l if k - min2 < 2.5]))
        for x, y in zip(target_luao_a, target_luao_b):
            if luaoa_l.count(x) > 1:
                at_x = -1
                for c in range(0, luaoa_l.count(x)):
                    at_x = luaoa_l.index(x, at_x + 1)
                    luaoa_d[i][at_x] = x
                    at_l3.append(at_x)
            else:
                at_x = luaoa_l.index(x)
                luaoa_d[i][at_x] = x
                at_l3.append(at_x)
            if luaob_l.count(y) > 1:
                at_y = -1
                for c in range(0, luaob_l.count(y)):
                    at_y = luaob_l.index(y, at_y + 1)
                    luaob_d[i][at_y] = x
                    at_l4.append(at_y)
            else:
                at_y = luaob_l.index(y)
                luaob_d[i][at_y] = y
                at_l4.append(at_y)

            at_l3 = list(set(at_l3))
            at_l4 = list(set(at_l4))

    hoaoa_df = pd.DataFrame(hoaoa_d, index=[i for i in at_l1], columns=[i for i in range(len(reactants))])
    hoaob_df = pd.DataFrame(hoaob_d, index=[i for i in at_l2], columns=[i for i in range(len(reactants))])
    luaoa_df = pd.DataFrame(luaoa_d, index=[i for i in at_l3], columns=[i for i in range(len(reactants))])
    luaob_df = pd.DataFrame(luaob_d, index=[i for i in at_l4], columns=[i for i in range(len(reactants))])

    # replace the missing data
    for i in [luaoa_df, luaob_df]:
        for j in i.columns:
            cl = list(i.loc[:, j])
            trigger = False
            for k in range(len(cl)):
                if type(cl[k]) == str:
                    i.iloc[0, j] = float('inf')
                else:
                    if not math.isnan(cl[k]):
                        trigger = True
            if not trigger and type(i.iloc[0, j]) != float:
                i.iloc[0, j] = -40
                print("fill missing data")
    for i in [hoaoa_df, hoaob_df]:
        for j in i.columns:
            cl = list(i.loc[:, j])
            trigger = False
            for k in range(len(cl)):
                if type(cl[k]) == str:
                    i.iloc[0, j] = float('inf')
                else:
                    if not math.isnan(cl[k]):
                        trigger = True
            if not trigger and type(i.iloc[0, j]) != float:
                i.iloc[0, j] = 40
                print("fill missing data")

    react_ss = {}
    react_ss = build_gap_d(react_ss, hoaoa_df, luaoa_df)
    react_ss = litgap(react_ss, hoaoa_df, luaoa_df)

    react_ss = build_gap_d(react_ss, hoaob_df, luaob_df)
    react_ss = litgap(react_ss, hoaob_df, luaob_df)


    react_os = {}
    react_os = build_gap_d(react_os, hoaoa_df, luaob_df)
    react_os = litgap(react_os, hoaoa_df, luaob_df)
    react_os = build_gap_d(react_os, hoaob_df, luaoa_df)
    react_os = litgap(react_os, hoaob_df, luaoa_df)
    return react_os, react_ss


def spr_intra_new(reactants, rec, hydrogen=True):
    dis_d = {}
    # use bo as measurement
    for smi in reactants:
        try:
            wboe = rec[smi]['wboe']
            mlfeatures = rec[smi]['ml_feature']
            if isinstance(wboe, pd.DataFrame):
            # wboe = pd.read_csv(f'{wd_pre}/{smi}/ml/wboe', header=None, sep='\s+')
            # if os.path.getsize(f'{wd_pre}/{smi}/ml/wboe') != 0:
            #     mlfeatures = pd.read_csv(f'{wd_pre}/{smi}/ml/ml_feature.csv')
                tmp = {}
                for i in range(wboe.shape[0]):
                    dagp = wboe.iloc[i, 3]
                    at_1, at_2 = wboe.iloc[i, 0], wboe.iloc[i, 1]
                    # try:
                    #     cn_1 = int(mlfeatures.loc[at_1 - 1, 'coordination number'])
                    # except:
                    #     cn_1 = 0
                    # try:
                    #     cn_2 = int(mlfeatures.loc[at_2 - 1, 'coordination number'])
                    # except:
                    #     cn_2 = 0
                    ele_1 = mlfeatures.loc[at_1 - 1, 'Atom']
                    ele_2 = mlfeatures.loc[at_2 - 1, 'Atom']
                    if not hydrogen and (int(ele_1) == 1 or int(ele_2) == 1):
                        pass
                    else:
                        tmp[dagp] = [smi, at_1, at_2]
                if not tmp:  # only h-bonds in struc
                    for i in range(wboe.shape[0]):
                        dagp = wboe.iloc[i, 3]
                        at_1, at_2 = wboe.iloc[i, 0], wboe.iloc[i, 1]

                        tmp[dagp] = [smi, at_1, at_2]
                dis_d.update(tmp)
        except:
            pass
            traceback.print_exc()
            # print(f'{wd_pre}/{smi}')
    return dis_d
