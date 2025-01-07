"""
spr = search possiable reaction site
@Author  ：Ying Chen
"""
import sys
import traceback
from itertools import product
import numpy as np
import pandas as pd
from aRST.script.dictionaryInspector import DictionaryInspector


def build_react_l(h_df, l_df, search_mode):
    """
    loop over Hdf and Ldf to build up atomic pairs
    """
    react_l = []
    for il, (istruc_l, row_l) in enumerate(l_df.items()):
        for ih, (istruc_h, row_h) in enumerate(h_df.items()):
            # Test search_mode
            if search_mode == 0 and il != ih:  # only search struc1 == struc2
                continue
            elif search_mode == 1 and il == ih:  # only search struc1 != struc2
                continue
    # for istruc_l, row_l in l_df.items():
    #     for istruc_h, row_h in h_df.items():
    #         # test search_mode
    #         if search_mode==0 and istruc_l!=istruc_h: # only search struc1 == struc2
    #             continue
    #         elif search_mode==1 and istruc_l==istruc_h: # only search struc1 != struc2
    #             continue

            # loop row to get H and L values
            for iat_l,L in  row_l.items():
                if isinstance(L, str) or np.isnan(L) or np.isinf(L):
                    continue
                for iat_h,H in  row_h.items():
                    if isinstance(H, str) or np.isnan(H) or np.isinf(H):
                        continue
                    gap = np.round(np.round(L, 2) - np.round(H, 2), 2)
                    if int(istruc_l)>=int(istruc_h):
                        react_l.append({"struc_combi":[str(istruc_h),str(istruc_l)],
                                        "at_combi":[str(iat_h),str(iat_l)],
                                        "HOAO":H,
                                        "LUAO": L,
                                        "gap": gap})
                    else:
                        react_l.append({"struc_combi": [str(istruc_l), str(istruc_h)],
                                        "at_combi": [str(iat_l), str(iat_h)],
                                        "HOAO": H,
                                        "LUAO": L,
                                        "gap": gap})
    return react_l


def refine_react_l(react_l,identical_atom_d):
    # merge degenrate pairs and only keep lowest gap
    out_react_l = []
    used_combi = []
    for item in react_l:
        struc1,struc2 = item["struc_combi"]
        at1,at2  = item["at_combi"]
        if tuple([struc1,struc2,(at1,at2)]) not in used_combi:
            # add all degenrated pair
            degenerate_at1=[at1]
            degenerate_at2 = [at2]
            for tmp in identical_atom_d[struc1]:
                if at1 in tmp["atom"]:
                    degenerate_at1 = tmp["atom"]
                    break
            for tmp in identical_atom_d[struc2]:
                if at2 in tmp["atom"]:
                    degenerate_at2 = tmp["atom"]
                    break
            all_combi = list(product(degenerate_at1, degenerate_at2))
            for combi in all_combi:
                used_combi.append(tuple([struc1,struc2,combi]))
            item["at_combi"]=all_combi
            out_react_l.append(item)
        else:
            #change gap
            for searchitem in out_react_l:
                search_allcombi = searchitem["at_combi"]
                if tuple([at1,at2]) in search_allcombi:
                    currentgap = item["gap"]
                    oldgap = searchitem["gap"]
                    if currentgap<oldgap:
                        searchitem["gap"]=currentgap
                        searchitem["HOAO"]=item["HOAO"]
                        searchitem["LUAO"] = item["LUAO"]
                        break
    return out_react_l



def SPR_HL(reactants,  # reactants id list
           struc_rec,  # struc dict record
           allow_hydrogen=True,  # if search range including hydrogen elements
           allow_halogen=True,  # if search range including halogen elements
           search_mode=1,  # 0 : only search intramollecule reaction possibility;
           # 1: only intermolecular;
           # 2: both inter- and intra reaction possibility.
           spr_mode=0 #0(default): globally compre both same spin and oppo spin HL combinations
                                    # 1: only compre same spin (HL_aa and HL_bb)
                                    # 2: only compre opposite spin (HL_ab and HL_ba)
           ):

    struc_inspector = DictionaryInspector(struc_rec)
    # collecting HL values into dataframe
    maxnat = max(struc_inspector.get_value(f'{i}.struc_info.reactivity.ml').shape[0] for i in reactants)
    HOAO_a = pd.DataFrame(columns=reactants, index=range(maxnat))
    HOAO_b = pd.DataFrame(columns=reactants, index=range(maxnat))
    LUAO_a = pd.DataFrame(columns=reactants, index=range(maxnat))
    LUAO_b = pd.DataFrame(columns=reactants, index=range(maxnat))
    newcolumn_names = ["HOAO_a (eV)", "HOAO_b (eV)", "LUAO_a (eV)", "LUAO_b (eV)"]
    dfs = [HOAO_a, HOAO_b, LUAO_a, LUAO_b]

    for istruc in range(len(reactants)):
        strucid = reactants[istruc]
        mldf_original = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.ml').copy()

        # remove all hydrogen/halogen if not allowed
        if not allow_hydrogen:
            # get index
            hindex = mldf_original[mldf_original['Atom'] == 1].index.tolist()
            mldf_original.loc[hindex, "HOAO_a (eV)"] = -np.inf
            mldf_original.loc[hindex, "HOAO_b (eV)"] = -np.inf
            mldf_original.loc[hindex, "LUAO_a (eV)"] = np.inf
            mldf_original.loc[hindex, "LUAO_a (eV)"] = np.inf

        if not allow_halogen:
            # get index
            hindex = mldf_original[mldf_original['Atom'] == 9 or
                                   mldf_original['Atom'] == 17 or
                                   mldf_original['Atom'] == 35 or
                                   mldf_original['Atom'] == 53 or
                                   mldf_original['Atom'] == 85].index.tolist()
            # change values
            mldf_original.loc[hindex, "HOAO_a (eV)"] = -np.inf
            mldf_original.loc[hindex, "HOAO_b (eV)"] = -np.inf
            mldf_original.loc[hindex, "LUAO_a (eV)"] = np.inf
            mldf_original.loc[hindex, "LUAO_a (eV)"] = np.inf

        for df, column_name in zip(dfs, newcolumn_names):
            tmpl = list(mldf_original.loc[:, column_name])
            for j in range(maxnat - len(tmpl)):
                tmpl.append(np.nan)
            df.iloc[:, istruc] = tmpl  # .reset_index(drop=True)

    identical_atom_d = find_identical_at(reactants,struc_rec)


    final_react = []
    if spr_mode==1:
        react_aa = build_react_l(HOAO_a, LUAO_a, search_mode)
        react_bb = build_react_l(HOAO_b, LUAO_b, search_mode)
        for react in [react_aa, react_bb]:
            final_react.extend(react)
    elif spr_mode==2:
        react_ab = build_react_l(HOAO_a, LUAO_b, search_mode)
        react_ba = build_react_l(HOAO_b, LUAO_a, search_mode)
        for react in [react_ab, react_ba]:
            final_react.extend(react)

    else:
        react_aa = build_react_l(HOAO_a, LUAO_a, search_mode)
        react_bb = build_react_l(HOAO_b, LUAO_b, search_mode)
        react_ab = build_react_l(HOAO_a, LUAO_b, search_mode)
        react_ba = build_react_l(HOAO_b, LUAO_a, search_mode)

        for react in [react_aa,react_bb,react_ab,react_ba]:
            final_react.extend(react)

    final_react = sorted(refine_react_l(final_react,identical_atom_d), key=lambda x: x['gap'])
    return final_react


def find_identical_at(reactants, struc_rec, tol=0.01):
    struc_inspector = DictionaryInspector(struc_rec)
    out = {}
    for istruc in range(len(reactants)):
        out[str(istruc)] = []
        strucid = reactants[istruc]
        mldf_original = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.ml').copy()
        hoao = mldf_original.loc[:, "HOAO_a (eV)"].tolist()
        hoaos_with_indices = sorted(enumerate(hoao), key=lambda x: x[1])
        indexed_groups = []
        current_group = [hoaos_with_indices[0]]
        for num_with_index in hoaos_with_indices[1:]:  # Iterate through the rest
            if abs(num_with_index[1] - current_group[-1][1]) <= tol:  # Check tolerance
                current_group.append(num_with_index)  # Add to current group
            else:
                indexed_groups.append(current_group)  # Save current group
                current_group = [num_with_index]  # Start a new group

        indexed_groups.append(current_group)  # Add the last group
        out[str(istruc)] = [{"atom": [str(item[0]) for item in group], "hoao": [item[1] for item in group]}
                       for group in indexed_groups]
    return out


def SPR_dissociation(reactants,  # reactants id list
                     struc_rec,  # struc dict record
                     allow_hydrogen=True,  # if search range including hydrogen elements
                     allow_halogen=True,  # if search range including halogen elements
                     search_mode=0  # 0: use wboe; 1: use MBO
                     ):
    struc_inspector = DictionaryInspector(struc_rec)
    BO = {}
    for strucid in reactants:
        wboe = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.bo')
        ml = struc_inspector.get_value(f'{strucid}.struc_info.reactivity.ml')
        if wboe is not None:
            for index in range(wboe.shape[0]):
                if search_mode == 1:
                    dgap = wboe.iloc[index, 3]  # MBO
                    if dgap <= 0.3:
                        continue
                else:
                    dgap = wboe.iloc[index, 2]  # energy-baesd MBO
                at1, at2 = wboe.iloc[index, :2]  # index starting from 1
                ele1, ele2 = ml.loc[at1 - 1, 'Atom'], ml.loc[at2 - 1, 'Atom']

                if not allow_hydrogen:
                    if int(ele1) == 1 or int(ele2) == 1:
                        continue
                if not allow_halogen:
                    if int(ele1) in [9, 17, 35, 53, 85] or int(ele2) in [9, 17, 35, 53, 85]:
                        continue
                value = [strucid, at1 - 1, at2 - 1, dgap]
                if dgap not in BO.keys():
                    BO[dgap] = [value]
                else:
                    BO[dgap].append(value)
    return BO


def build_bo_d(reactants,
               struc_rec,  # struc dict record
               allow_hydrogen=True,  # if search range including hydrogen elements
               allow_halogen=True,  # if search range including halogen elements
               ):
    struc_inspector = DictionaryInspector(struc_rec)
    bo = {}
    for strucid in reactants:
        bo[strucid] = []
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
                    if not allow_halogen:
                        if int(ele_1) in [9, 17, 35, 53, 85] or int(ele_2) in [9, 17, 35, 53, 85]:
                            continue
                    value = [f'{at_1 - 1}', f'{at_2 - 1}']
                    bo[strucid].append(value)
    return bo

