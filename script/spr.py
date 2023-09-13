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


# def litgap(react_d, hoao_df, luao_df):
#
#     tmp1 = luao_df.min().tolist()
#     tmp2 = hoao_df.max().tolist()
#     tmp2c = copy.deepcopy(tmp2)
#     tmp2c.pop(tmp1.index(min(tmp1)))
#     lit1 = min(tmp1) - max(tmp2c)
#
#     tmp1c = copy.deepcopy(tmp1)
#     tmp1c.pop(tmp2.index(max(tmp2)))
#     lit2 = min(tmp1c) - max(tmp2c)
#     del_l = [i for i in react_d.keys() if i < min(lit1, lit2)]
#     for i in del_l:
#         del react_d[i]
#     return react_d


def fill_miss_data(l, typ):
    for i in l:
        for j in i:
            if type(j) == str or math.isnan(j):
                if typ == 'h':
                    i[i.index(j)] = -float('inf')
                else:
                    i[i.index(j)] = float('inf')

    return l


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

    # merge equivalent pairs
    result_dict = {}
    for gap in react_d.keys():
        value = react_d[gap]
        # print(f'gap = {gap}, value = {value}')
        gap_appoxi = round(gap,4)
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
                if at1 not in item and at2 not in item:
                    if [at1, at2] not in at_l and [at2, at1] not in at_l:
                        if struc1 == struc1_old and struc2 == struc2_old:
                            at_l.append([at1, at2])
                            value_old[1] = str(at_l)
                            result_dict[gap_appoxi] = value_old
                        elif struc2 == struc1_old and struc1 == struc2_old:
                            at_l.append([at2, at1])
                            value_old[1] = str(at_l)
                            result_dict[gap_appoxi] = value_old

#     print("hier",result_dict)
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

    # find max in hoao and min in luao
    if 'hoao_a' in list(rec[reactants[0]]['ml_feature'].columns):
        HOAO_a = [list(rec[k]['ml_feature']['hoao_a']) for k in reactants]
        HOAO_b = [list(rec[k]['ml_feature']['hoao_b']) for k in reactants]
        LUAO_a = [list(rec[k]['ml_feature']['luao_a']) for k in reactants]
        LUAO_b = [list(rec[k]['ml_feature']['luao_b']) for k in reactants]
    else:
        HOAO_a = [list(rec[k]['ml_feature']['hoao']) for k in reactants]
        HOAO_b = [list(rec[k]['ml_feature']['hoao']) for k in reactants]
        LUAO_a = [list(rec[k]['ml_feature']['luao']) for k in reactants]
        LUAO_b = [list(rec[k]['ml_feature']['luao']) for k in reactants]
    # replace the missing data
    HOAO_a = fill_miss_data(HOAO_a, 'h')
    HOAO_b = fill_miss_data(HOAO_b, 'h')
    LUAO_a = fill_miss_data(LUAO_a, 'l')
    LUAO_b = fill_miss_data(LUAO_b, 'l')

    max_a = max([max(k) for k in HOAO_a])
    max_b = max([max(k) for k in HOAO_b])
    min_a = min([min(k) for k in LUAO_a])
    min_b = min([min(k) for k in LUAO_b])

    # form target atom df
    for i in range(len(reactants)):
        hoaoa_d[i] = {}
        hoaob_d[i] = {}
        luaoa_d[i] = {}
        luaob_d[i] = {}

        hoaoa_l = HOAO_a[i]
        hoaob_l = HOAO_b[i]
        luaoa_l = LUAO_a[i]
        luaob_l = LUAO_b[i]

        lim = 5
        target_hoao_a = [k for k in hoaoa_l if max_a - k < lim]
        target_hoao_b = [k for k in hoaob_l if max_b - k < lim]
        target_luao_a = [k for k in luaoa_l if k - min_a < lim]
        target_luao_b = [k for k in luaob_l if k - min_b < lim]
        for x in target_hoao_a:
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

        for y in target_hoao_b:
            if hoaob_l.count(y) > 1:
                at_y = -1
                for c in range(0, hoaob_l.count(y)):
                    at_y = hoaob_l.index(y, at_y + 1)
                    hoaob_d[i][at_y] = y
                    at_l2.append(at_y)
            else:
                at_y = hoaob_l.index(y)
                hoaob_d[i][at_y] = y
                at_l2.append(at_y)
        at_l1 = list(set(at_l1))
        at_l2 = list(set(at_l2))

        for x in target_luao_a:
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
        for y in target_luao_b:
            if luaob_l.count(y) > 1:
                at_y = -1
                for c in range(0, luaob_l.count(y)):
                    at_y = luaob_l.index(y, at_y + 1)
                    luaob_d[i][at_y] = y
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

    print(hoaoa_df)
    print(hoaob_df)
    print(luaoa_df)
    print(luaob_df)
    react_ss = {}
    react_ss = build_gap_d(react_ss, hoaoa_df, luaoa_df)
    # react_ss = litgap(react_ss, hoaoa_df, luaoa_df)

    react_ss = build_gap_d(react_ss, hoaob_df, luaob_df)

    # react_ss = litgap(react_ss, hoaob_df, luaob_df)


    react_os = {}
    react_os = build_gap_d(react_os, hoaoa_df, luaob_df)
    # react_os = litgap(react_os, hoaoa_df, luaob_df)
    react_os = build_gap_d(react_os, hoaob_df, luaoa_df)
    # react_os = litgap(react_os, hoaob_df, luaoa_df)

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
