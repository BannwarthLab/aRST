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
        LUAO_a = [list(rec[k]['ml_feature']['LUAO_a (eV)']) for k in reactants]
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

        lim = 2.5
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