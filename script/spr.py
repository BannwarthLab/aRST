import traceback
import numpy as np
import pandas as pd


def build_gap_d(react_d, h_df, l_df):
    # loop 2 df to get HL gap
    for index1, row1 in l_df.iterrows():
        for column1 in l_df.columns:
            if pd.notna(l_df.at[index1, column1]):
                L = l_df.at[index1, column1]
                for index2, row2 in h_df.iterrows():
                    for column2 in h_df.columns:
                        if pd.notna(h_df.at[index2, column2]):
                            H = h_df.at[index2, column2]
                            if isinstance(H, float) and isinstance(L, float):
                                gap = np.round(L - H, 3)
                                struc1, struc2 = str(column2), str(column1)
                                tmp_key = (struc1, struc2)
                                tmp_val = [str(index2), str(index1)]
                                # check if gap exist
                                if gap not in react_d.keys():
                                    react_d[gap] = {}
                                    react_d[gap][tmp_key] = [tmp_val]
                                # degenrated pairs form different struc combinations
                                elif tmp_key not in react_d[gap].keys() \
                                        and tuple(list(tmp_key)[::-1]) not in react_d[gap].keys():
                                    react_d[gap][tmp_key] = [tmp_val]
                                # degenrated pairs form same struc combinations
                                else:
                                    for current_key, current_item in react_d[gap].items():
                                        if tmp_key == current_key and len(set(tmp_key)) == 1:
                                            if tmp_val in current_item or tmp_val[::-1] in current_item:
                                                continue
                                        if tmp_key == current_key and tmp_val \
                                                not in current_item:
                                            # same atom not allow to appear twice e.g. [['2', '1'], ['3', '1']]
                                            if tmp_val[0] not in [k[0] for k in current_item] \
                                                    and tmp_val[1] not in [k[1] for k in current_item]:
                                                current_item.append(tmp_val)
                                        elif tuple(list(tmp_key)[::-1]) == current_key and tmp_val[::-1] \
                                                not in current_item:
                                            if tmp_val[0] not in [k[0] for k in current_item] \
                                                    and tmp_val[1] not in [k[1] for k in current_item]:
                                                current_item.append(tmp_val[::-1])

                            ## text-wise out put
                            # value = f"1|{column2},{column1}|{index2},{index1}"
                            #
                            # # check if gap exist
                            # if gap not in react_d:
                            #     react_d[gap] = value
                            # else:
                            #     # check if gap degenerated
                            #     old_value = react_d[gap].split("|")
                            #     for struc_combi in [old_value[i] for i in range(1, len(old_value), 2)]:
                            #         a, b = struc_combi.split(",")
                            #         str_atl = old_value[old_value.index(struc_combi) + 1]
                            #         at_l = str_atl.split(",")
                            #         sep_atl = [at_l[k:k + 2] for k in range(0, len(at_l), 2)]
                            #         lpart = old_value[1:old_value.index(struc_combi) + 1]
                            #         rpart = old_value[old_value.index(struc_combi) + 2:]
                            #         tmp = None
                            #         if int(column2) == int(a) and int(column1) == int(b):
                            #             if [f"{index2}", f"{index1}"] not in sep_atl:
                            #                 tmp = f"{index2},{index1}"
                            #         elif int(column2) == int(b) and int(column1) == int(a):
                            #             if [f"{index1}", f"{index2}"] not in sep_atl:
                            #                 tmp = f"{index1},{index2}"
                            #         if tmp is not None:
                            #             new_count = int(old_value[0]) + 1
                            #             new_value = [str(new_count)] + lpart + [str_atl + "," + tmp] + rpart
                            #             new_value = "|".join(new_value)
                            #             react_d[gap] = new_value
    return react_d


def SPR_HL(reactants, struc_inspector, allow_hydrogen=True):
    """
    :param allow_hydrogen:
    :param reactants: strucid of reactants
    :param struc_inspector: struc_record
    :return: d[gap] = "degenerate number (int)|strucid1(without 0),strucid2|at1,at2|....."
    """
    # Calculate the maximum number of rows among Reactant molecules
    nat = max(struc_inspector.get_value(f'{i}.struc_info.reactivity.ml').shape[0] for i in reactants)

    # Initialize DataFrames with NaN values
    HOAO_a = pd.DataFrame(index=range(nat))
    HOAO_b = pd.DataFrame(index=range(nat))
    LUAO_a = pd.DataFrame(index=range(nat))
    LUAO_b = pd.DataFrame(index=range(nat))

    # Define column names
    column_names = ["HOAO_a (eV)", "HOAO_b (eV)", "LUAO_a (eV)", "LUAO_b (eV)"]
    dfs = [HOAO_a, HOAO_b, LUAO_a, LUAO_b]

    # build HOAO_df and LUAO_df
    if 'HOAO_a (eV)' in struc_inspector.get_value(f'{reactants[0]}.struc_info.reactivity.ml').columns:
        for df, column_name in zip(dfs, column_names):
            for strucid in reactants:
                tmpdf = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.ml').copy(deep=True)
                # remove all hydrogen if not allowed
                if not allow_hydrogen:
                    for at in range(tmpdf.shape[0]):
                        if int(tmpdf.iloc[at]['Atom']) == 1:
                            tmpdf.iloc[at][f"{column_name}"] = None
                df[strucid.lstrip("0")] = tmpdf[column_name].reset_index(drop=True)
            # Convert all elements in the DataFrame to numeric, coercing errors to NaN
            df.apply(pd.to_numeric, errors='coerce')

    # remove all Value exceeding the limit
    lim = 200  # eV

    for df in [HOAO_a, HOAO_b]:
        max_value = pd.to_numeric(df.values.flatten(), errors='coerce').max()
        # df[df < max_value - lim] = None
        df.mask((df < max_value - lim) | df.applymap(lambda x: isinstance(x, str)), inplace=True)

    for df in [LUAO_a, LUAO_b]:
        min_value = pd.to_numeric(df.values.flatten(), errors='coerce').min()
        df.mask((df > min_value + lim) | df.applymap(lambda x: isinstance(x, str)), inplace=True)
        # df[df > min_value + lim] = None


    # buid HL_d
    react_ss = {}
    react_ss = build_gap_d(react_ss, HOAO_a, LUAO_a)
    react_ss = build_gap_d(react_ss, HOAO_b, LUAO_b)

    react_os = {}
    react_os = build_gap_d(react_os, HOAO_a, LUAO_b)
    react_os = build_gap_d(react_os, HOAO_b, LUAO_a)

    # for print out
    HOAO_a.index.name = 'HOAO_a'
    HOAO_b.index.name = 'HOAO_b'
    LUAO_a.index.name = 'LUAO_a'
    LUAO_b.index.name = 'LUAO_b'
    return react_os, react_ss, HOAO_a, HOAO_b, LUAO_a, LUAO_b


def SPR_dissociation(reactants, struc_inspector, allow_hydrogen=True,allow_halogen=True,type="bo"):
    """
    :param reactants: strucid of reactants
    :param struc_inspector: struc_record
    :param allow_hydrogen: count H-bonds in or not
    :return: d[bo]:
    """
    BO = {}
    for strucid in reactants:
        wboe = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.bo')
        ml = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.ml')
        if wboe is not None:
            for index in range(wboe.shape[0]):
                if type == "bo":
                    dgap = wboe.iloc[index, 3]
                else:
                    dgap = wboe.iloc[index, 2]
                if 0.3 < dgap if type == "bo" else dgap:
                    at_1, at_2 = wboe.iloc[index, 0], wboe.iloc[index, 1]
                    ele_1, ele_2 = ml.iloc[at_1-1]['Atom'], ml.iloc[at_2-1]['Atom']
                    if not allow_hydrogen and (int(ele_1) == 1 or int(ele_2) == 1):
                        continue
                    if not allow_halogen and (int(ele_1) in [9,17,35,53,85] or int(ele_2) in [9,17,35,53,85]):
                        continue
                    value = [[strucid.lstrip("0"),str(at_1 - 1), str(at_2 - 1), wboe.iloc[index, 2]]]
                    if dgap not in BO.keys():
                        BO[dgap] = value
                    else:
                        BO[dgap].append([strucid.lstrip("0"), str(at_1 - 1), str(at_2 - 1), wboe.iloc[index, 2]])

                    # value = f'{strucid.lstrip("0")}|{at_1 - 1},{at_2 - 1}|{wboe.iloc[index, 2]}'
                    # if dgap not in BO.keys():
                    #     BO[dgap] = value
                    # else:
                    #     oldvalue = BO[dgap]
                    #     BO[dgap] = oldvalue + "|" + value

    return BO


# for react_d in [react_os, react_ss]:
#     for current_key, current_item in react_d.items():
#         # 5.93, 1|1,1|0,0
#         item_l = current_item.split("|")
#         if len(set(item_l[1].split(","))) == 1:  # intramolecular
#             struc = item_l[1].split(",")[0]
#             ndegenerate = int(item_l[0])
#             str_atl = item_l[2]  # only consider first molecule?
#             at_l = str_atl.split(",")
#             sep_atl = [at_l[k:k + 2] for k in range(0, len(at_l), 2)]
#             if key is None or current_key < key:  # low HL_gap
#                 key = current_key
#                 push_struc = struc
#                 push_atl = sep_atl


def build_bo_d(reactants, struc_inspector, allow_hydrogen=True):
    bo = {}
    for strucid in reactants:
        bo[strucid.lstrip("0")] = []
        wboe = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.bo')
        if wboe is not None:
            ml = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.ml')
            for index in range(wboe.shape[0]):
                dgap = wboe.iloc[index, 3]
                if 0.3 < dgap:
                    at_1, at_2 = wboe.iloc[index, 0], wboe.iloc[index, 1]
                    ele_1, ele_2 = ml.iloc[at_1 - 1]['Atom'], ml.iloc[at_2 - 1]['Atom']
                    if not allow_hydrogen and (int(ele_1) == 1 or int(ele_2) == 1):
                        continue
                    value = [f'{at_1-1}', f'{at_2-1}']
                    bo[strucid.lstrip("0")].append(value)
    return bo


def spr_intra_new(reactants, rec, mode=0, hydrogen=True):
    dis_d = {}
    if mode == 0:  # use bo as measurement
        colunmnr = 3
    else:
        colunmnr = 2  # use bde as measurement

    for smi in reactants:
        try:
            wboe = rec[smi]['wboe']
            mlfeatures = rec[smi]['ml_feature']
            if isinstance(wboe, pd.DataFrame):
                tmp = {}
                for i in range(wboe.shape[0]):
                    dagp = wboe.iloc[i, colunmnr]
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
                    elif colunmnr == 3 and (dagp <= 0.3 or dagp > 1.2):
                        pass
                    else:
                        tmp[dagp] = [smi, at_1, at_2]
                # if not tmp:  # only h-bonds in struc
                #     for i in range(wboe.shape[0]):
                #         dagp = wboe.iloc[i, colunmnr]
                #         at_1, at_2 = wboe.iloc[i, 0], wboe.iloc[i, 1]
                #         tmp[dagp] = [smi, at_1, at_2]
                dis_d.update(tmp)
        except:
            pass
            traceback.print_exc()
            # print(f'{wd_pre}/{smi}')
    return dis_d


# def spr_inter_new(reactants, rec):
#     hoaoa_d = {}
#     at_l1 = []
#     hoaob_d = {}
#     at_l2 = []
#
#     luaoa_d = {}
#     at_l3 = []
#     luaob_d = {}
#     at_l4 = []
#
#     # find max in hoao and min in luao
#     if 'HOAO_a (eV)' in list(rec[reactants[0]]['ml_feature'].columns):
#         HOAO_a = [list(rec[k]['ml_feature']['HOAO_a (eV)']) for k in reactants]
#         HOAO_b = [list(rec[k]['ml_feature']['HOAO_b (eV)']) for k in reactants]
#         LUAO_a = [list(rec[k]['ml_feature']['LUAO_a (eV)']) for k in reactants]
#         LUAO_b = [list(rec[k]['ml_feature']['LUAO_b (eV)']) for k in reactants]
#     else:
#         HOAO_a = [list(rec[k]['ml_feature']['HOAO (eV)']) for k in reactants]
#         HOAO_b = [list(rec[k]['ml_feature']['HOAO (eV)']) for k in reactants]
#         LUAO_a = [list(rec[k]['ml_feature']['LUAO (eV)']) for k in reactants]
#         LUAO_b = [list(rec[k]['ml_feature']['LUAO (eV)']) for k in reactants]
#     # replace the missing data
#     HOAO_a = fill_miss_data(HOAO_a, 'h')
#     HOAO_b = fill_miss_data(HOAO_b, 'h')
#     LUAO_a = fill_miss_data(LUAO_a, 'l')
#     LUAO_b = fill_miss_data(LUAO_b, 'l')
#
#     max_a = max([max(k) for k in HOAO_a])
#     max_b = max([max(k) for k in HOAO_b])
#     min_a = min([min(k) for k in LUAO_a])
#     min_b = min([min(k) for k in LUAO_b])
#
#     # form target atom df
#     for i in range(len(reactants)):
#         hoaoa_d[i] = {}
#         hoaob_d[i] = {}
#         luaoa_d[i] = {}
#         luaob_d[i] = {}
#
#         hoaoa_l = HOAO_a[i]
#         hoaob_l = HOAO_b[i]
#         luaoa_l = LUAO_a[i]
#         luaob_l = LUAO_b[i]
#
#         lim = 500
#         # lim = 5
#         target_hoao_a = [k for k in hoaoa_l if max_a - k < lim]
#         target_hoao_b = [k for k in hoaob_l if max_b - k < lim]
#         target_luao_a = [k for k in luaoa_l if k - min_a < lim]
#         target_luao_b = [k for k in luaob_l if k - min_b < lim]
#         for x in target_hoao_a:
#             if hoaoa_l.count(x) > 1:
#                 at_x = -1
#                 for c in range(0, hoaoa_l.count(x)):
#                     at_x = hoaoa_l.index(x, at_x + 1)
#                     hoaoa_d[i][at_x] = x
#                     at_l1.append(at_x)
#             else:
#                 at_x = hoaoa_l.index(x)
#                 hoaoa_d[i][at_x] = x
#                 at_l1.append(at_x)
#
#         for y in target_hoao_b:
#             if hoaob_l.count(y) > 1:
#                 at_y = -1
#                 for c in range(0, hoaob_l.count(y)):
#                     at_y = hoaob_l.index(y, at_y + 1)
#                     hoaob_d[i][at_y] = y
#                     at_l2.append(at_y)
#             else:
#                 at_y = hoaob_l.index(y)
#                 hoaob_d[i][at_y] = y
#                 at_l2.append(at_y)
#         at_l1 = list(set(at_l1))
#         at_l2 = list(set(at_l2))
#
#         for x in target_luao_a:
#             if luaoa_l.count(x) > 1:
#                 at_x = -1
#                 for c in range(0, luaoa_l.count(x)):
#                     at_x = luaoa_l.index(x, at_x + 1)
#                     luaoa_d[i][at_x] = x
#                     at_l3.append(at_x)
#             else:
#                 at_x = luaoa_l.index(x)
#                 luaoa_d[i][at_x] = x
#                 at_l3.append(at_x)
#         for y in target_luao_b:
#             if luaob_l.count(y) > 1:
#                 at_y = -1
#                 for c in range(0, luaob_l.count(y)):
#                     at_y = luaob_l.index(y, at_y + 1)
#                     luaob_d[i][at_y] = y
#                     at_l4.append(at_y)
#             else:
#                 at_y = luaob_l.index(y)
#                 luaob_d[i][at_y] = y
#                 at_l4.append(at_y)
#
#         at_l3 = list(set(at_l3))
#         at_l4 = list(set(at_l4))
#
#     hoaoa_df = pd.DataFrame(hoaoa_d, index=[i for i in at_l1], columns=[i for i in range(len(reactants))])
#     hoaob_df = pd.DataFrame(hoaob_d, index=[i for i in at_l2], columns=[i for i in range(len(reactants))])
#     luaoa_df = pd.DataFrame(luaoa_d, index=[i for i in at_l3], columns=[i for i in range(len(reactants))])
#     luaob_df = pd.DataFrame(luaob_d, index=[i for i in at_l4], columns=[i for i in range(len(reactants))])
#
#     react_ss = {}
#     react_ss = build_gap_d(react_ss, hoaoa_df, luaoa_df)
#     # react_ss = litgap(react_ss, hoaoa_df, luaoa_df)
#
#     react_ss = build_gap_d(react_ss, hoaob_df, luaob_df)
#
#     # react_ss = litgap(react_ss, hoaob_df, luaob_df)
#
#     react_os = {}
#     react_os = build_gap_d(react_os, hoaoa_df, luaob_df)
#     # react_os = litgap(react_os, hoaoa_df, luaob_df)
#     react_os = build_gap_d(react_os, hoaob_df, luaoa_df)
#     # react_os = litgap(react_os, hoaob_df, luaoa_df)
#
#     # for print out
#     for i in [hoaoa_df, hoaob_df, luaoa_df, luaob_df]:
#         i.columns = [rec[reactants[i]]['mb'].split('|')[2] for i in range(len(reactants))]
#     hoaoa_df.index.name = 'HOAO_a'
#     hoaob_df.index.name = 'HOAO_b'
#     luaoa_df.index.name = 'LUAO_a'
#     luaob_df.index.name = 'LUAO_b'
#     return react_os, react_ss, hoaoa_df, hoaob_df, luaoa_df, luaob_df


def spr_intraassoci(mlfeatures, wboe):
    bonds = [[wboe.iloc[i, 0], wboe.iloc[i, 1]] for i in range(wboe.shape[0]) if 0.6 < wboe.iloc[i, 3]]
    intraasso_os = {}
    intraasso_ss = {}
    ha = list(mlfeatures['HOAO_a (eV)'])
    hb = list(mlfeatures['HOAO_b (eV)'])
    la = list(mlfeatures['LUAO_a (eV)'])
    lb = list(mlfeatures['LUAO_b (eV)'])
    for i in range(len(ha)):
        for j in range(len(ha)):
            if i != j:
                # simple replace if there is degenerated pair
                intraasso_ss[np.round(la[i] - ha[j], 3)] = [i + 1, j + 1]
                intraasso_ss[np.round(lb[i] - hb[j], 3)] = [i + 1, j + 1]
                intraasso_os[np.round(la[i] - hb[j], 3)] = [i + 1, j + 1]
                intraasso_os[np.round(lb[i] - ha[j], 3)] = [i + 1, j + 1]

    tmpd1 = intraasso_os.copy()
    for i in list(tmpd1.keys()):
        value = tmpd1[i]
        if value in bonds or reversed(value) in bonds:
            del intraasso_os[i]
    tmpd2 = intraasso_ss.copy()
    for i in list(tmpd2.keys()):
        value = tmpd2[i]
        if value in bonds or reversed(value) in bonds:
            del intraasso_ss[i]
    return intraasso_os, intraasso_ss
