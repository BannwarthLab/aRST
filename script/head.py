import os
import traceback

import pandas as pd

from Element import element2


def print_header():
    print('==========================================================================================')
    print('|                                                                                        |')
    print('|                         Reaction Path Serach Tool                                      |')
    print('|                                                                                        |')
    print('|                                 Y.Chen                                                 |')
    print('|                                                                                        |')
    print('|                      Institute of Physical Chemistry                                   |')
    print('|                          RWTH Aachen University                                        |')
    print('|                                                                                        |')
    print('|                           version : 1.0_alpha                                          |')
    print('|                                                                                        |')
    print('==========================================================================================')
    print('\n')


def call_header(text, typ="*"):
    print(text.center(90, typ))


def printout_cycle_results(path_i, reactants, association_d,dissociation_d,
                           hoaoa_df, hoaob_df, luaoa_df, luaob_df, rec):
    try:
        call_header(f"REACTANTS INFO FOR PATH {path_i}: 1-{len(association_d.keys()) + len(dissociation_d.keys())} ")
        # print(f"Reactants for this path {path_i}: {reactants}\n")
        for i in reactants:
            name = rec[i]['name']
            at_df = {}
            at_l = rec[i]['ml_feature'][['Atom']].to_dict()
            for j in list(at_l['Atom'].keys()):
                old_v = at_l['Atom'][j]
                new_v = element2[old_v]
                at_l['Atom'][j] = new_v
            at_l['dir'] = i
            at_df[name] = at_l
            call_header("Element", " ")
            print(at_df, '\n')
            call_header("WBOE", " ")
            wboe = rec[i]['wboe']
            if not isinstance(wboe, list):
                wboe.columns = ['At1', 'AT2', 'BDE', 'BO']
                print(wboe.to_string(index=False), '\n')

        call_header("ATOMISTIC FEATURES", " ")

        # Reset the indices of the dataframes
        df1_reset = hoaoa_df.reset_index()
        df2_reset = hoaob_df.reset_index()
        merged_df = pd.concat([df1_reset.round(3), df2_reset.round(3)], axis=1, keys=['HOAO_a', 'HOAO_b'])
        print(merged_df.to_string(index=False), '\n')

        df1_reset = luaoa_df.reset_index()
        df2_reset = luaob_df.reset_index()
        merged_df = pd.concat([df1_reset.round(3), df2_reset.round(3)], axis=1, keys=['LUAO_a', 'LUAO_b'])
        print(merged_df.to_string(index=False), '\n')

        # print(hoaoa_df.round(3), '\n')
        # print(hoaob_df.round(3), '\n')
        # print(luaoa_df.round(3), '\n')
        # print(luaob_df.round(3), '\n')

        call_header("REACTION INFO"," ")
        print("association:")
        for i in list(association_d.keys()):
            print(f"{i}: {association_d[i]} ")
        print("dissociation:")
        for i in list(dissociation_d.keys()):
            print(f"{i}: {dissociation_d[i]} ")


        # for i in sorted(push_dict_1.keys())[:5]:
        #     print(f'dict1: gap: {i},at_l: {push_dict_1[i]}\n')
        # for i in sorted(push_dict_2.keys())[:5]:
        #     print(f'dict2: gap: {i},at_l: {push_dict_2[i]}\n')
    except:
        traceback.print_exc()
    return
