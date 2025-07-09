import sys
from itertools import product
import numpy as np
import pandas as pd
import aRST.script.initialization as storage


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

            # loop row to get H and L values
            for iat_l, L in row_l.items():
                if isinstance(L, str) or np.isnan(L) or np.isinf(L):
                    continue
                for iat_h, H in row_h.items():
                    if isinstance(H, str) or np.isnan(H) or np.isinf(H):
                        continue
                    gap = np.round(np.round(L, 2) - np.round(H, 2), 2)
                    if int(istruc_l) <= int(istruc_h):
                        react_l.append({"struc_combi": [str(istruc_h), str(istruc_l)],
                                        "at_combi": [str(iat_h), str(iat_l)],
                                        "HOAO": H,
                                        "LUAO": L,
                                        "gap": gap})
                    else:
                        react_l.append({"struc_combi": [str(istruc_l), str(istruc_h)],
                                        "at_combi": [str(iat_l), str(iat_h)],
                                        "HOAO": H,
                                        "LUAO": L,
                                        "gap": gap})
    return react_l



def find_identical_at_from_afo(HOAO_a, HOAO_b, LUAO_a, LUAO_b, tol=0.01):
    out = {}

    # Process each structure that exists in HOAO_a
    for struc in HOAO_a.columns.unique():
        out[str(struc)] = []
        used_indices = set()

        for idx in HOAO_a.index:
            if idx in used_indices:
                continue

            # Get reference values
            ref_values = {
                'ha': HOAO_a.loc[idx, struc],
                'hb': HOAO_b.loc[idx, struc],
                'la': LUAO_a.loc[idx, struc],
                'lb': LUAO_b.loc[idx, struc]
            }

            # Ensure scalar values (handle Series cases)
            ref_values = {k: v.iloc[0] if isinstance(v, pd.Series) else v
                          for k, v in ref_values.items()}

            # Find degenerate atoms
            degenerate = [idx]
            used_indices.add(idx)

            # Compare with subsequent atoms
            for i in range(idx + 1, len(HOAO_a)):
                if i in used_indices:
                    continue

                # Get comparison values
                curr_values = {
                    'ha': HOAO_a.loc[i, struc],
                    'hb': HOAO_b.loc[i, struc],
                    'la': LUAO_a.loc[i, struc],
                    'lb': LUAO_b.loc[i, struc]
                }
                curr_values = {k: v.iloc[0] if isinstance(v, pd.Series) else v
                               for k, v in curr_values.items()}

                # Check tolerance
                if all(abs(curr_values[k] - ref_values[k]) < tol for k in ref_values):
                    degenerate.append(i)
                    used_indices.add(i)

            out[str(struc)].append({"atom": degenerate})

    return out


def refine_react_l(react_l, identical_atom_d, intra=False):
    # merge identical pais and keep the lowest gap
    # in the degenearted pairs, at can appear only once
    out_react_l = []
    used_combi = []
    if not intra:
        for item in react_l:
            struc1, struc2 = item["struc_combi"]  # istruc1 > istruc2
            at1, at2 = item["at_combi"]
            if tuple((int(struc1), int(struc2), int(at1), int(at2))) not in used_combi:
                # find all degenerated pair
                degenerate_at1 = [at1]
                degenerate_at2 = [at2]
                for tmp in identical_atom_d[struc1]:
                    if int(at1) in tmp["atom"]:
                        degenerate_at1 = tmp["atom"]
                        break
                for tmp in identical_atom_d[struc2]:
                    if int(at2) in tmp["atom"]:
                        degenerate_at2 = tmp["atom"]
                        break

                all_combi = list(product(degenerate_at1, degenerate_at2))
                if int(struc1) ==int(struc2):
                    all_combi2 = list(product(degenerate_at2, degenerate_at1))
                    all_combi.extend(all_combi2)
                used_l1 = []
                used_l2 = []
                pairs = []
                for combi in all_combi:
                    tat1, tat2 = combi
                    used_combi.append(tuple((int(struc1), int(struc2), int(tat1), int(tat2))))
                    if tat1 not in used_l1 and tat2 not in used_l2:
                        pairs.append([tat1, tat2])
                        used_l1.append(tat1)
                        used_l2.append(tat2)
                        # gap_value = next((tmpitem["gap"] for tmpitem in react_l if
                        #                   tmpitem["struc_combi"] == [str(struc1), str(struc2)]
                        #                   and tmpitem["at_combi"] == [str(tat1),str(tat2)]),None)

                item["at_combi"] = pairs

                if item not in out_react_l:
                    out_react_l.append(item)
            else:
                # change gap
                for searchitem in out_react_l:
                    search_struc1, search_struc2 = searchitem["struc_combi"]
                    if int(search_struc1) ==int(struc1) and int(search_struc2) ==int(struc2):
                        search_allcombi = searchitem["at_combi"]
                        triger = next((True for tmpitem in search_allcombi
                                       if int(tmpitem[0]) ==int(at1) and int(tmpitem[1]) ==int(at2)),False)
                        if triger:
                            currentgap = item["gap"]
                            oldgap = searchitem["gap"]
                            if currentgap < oldgap:
                                searchitem["gap"] = currentgap
                                searchitem["HOAO"] = item["HOAO"]
                                searchitem["LUAO"] = item["LUAO"]
                                break
    else:
        # todo
        for item in react_l:
            struc = item["struc_combi"][0]
            at1, at2 = item["at_combi"]
            if tuple((int(struc), int(at1), int(at2))) not in used_combi:
                # find all degenerated pair
                degenerate_at1 = [at1]
                degenerate_at2 = [at2]
                for tmp in identical_atom_d[struc]:
                    if at1 in tmp["atom"]:
                        degenerate_at1 = tmp["atom"]
                        break
                for tmp in identical_atom_d[struc]:
                    if at2 in tmp["atom"]:
                        degenerate_at2 = tmp["atom"]
                        break
                all_combi = list(product(degenerate_at1, degenerate_at2))
                for combi in all_combi:
                    tat1, tat2 = combi
                    used_combi.append(tuple((int(struc), int(tat1), int(tat2))))
                used_l = []
                pairs = []
                for pair in all_combi:
                    at1, at2 = pair
                    if at1 not in used_l and at2 not in used_l:
                        pairs.append([at1, at2])
                        used_l.append(at1)
                        used_l.append(at2)
                item["at_combi"] = pairs
                if item not in out_react_l:
                    out_react_l.append(item)
            else:
                # change gap
                for searchitem in out_react_l:
                    search_struc = searchitem["struc_combi"][0]
                    if search_struc == struc:
                        search_allcombi = searchitem["at_combi"]
                        if tuple([at1, at2]) in search_allcombi:
                            currentgap = item["gap"]
                            oldgap = searchitem["gap"]
                            if currentgap < oldgap:
                                searchitem["gap"] = currentgap
                                searchitem["HOAO"] = item["HOAO"]
                                searchitem["LUAO"] = item["LUAO"]
                                break
        # change key name struc_combi to struc
        for item in out_react_l:
            item["struc"] = item["struc_combi"][0]
            del item["struc_combi"]
            # make sure iat1>iat2
            pairs = [[max(at1, at2), min(at1, at2)] for at1, at2 in item["at_combi"]]
            item["at_combi"] = pairs
    return out_react_l


class SPR:
    # atom index starts with 0
    def __init__(self, reactants,
                 allow_hydrogen=True,
                 allow_halogen=True,
                 afo_l=None):
        self.reactants = reactants
        self.allowhydrogen = allow_halogen
        self.allowhalogen = allow_halogen

        if afo_l is None:
            maxnat = max(storage.allstruc.get_value(f'{i}.struc_info.reactivity.afo').shape[0] for i in self.reactants)
        else:
            maxnat = max(i.shape[0] for i in afo_l)

        # collecting HL values into dataframe
        self.HOAO_a = pd.DataFrame(columns=self.reactants, index=range(maxnat))
        self.HOAO_b = pd.DataFrame(columns=self.reactants, index=range(maxnat))
        self.LUAO_a = pd.DataFrame(columns=self.reactants, index=range(maxnat))
        self.LUAO_b = pd.DataFrame(columns=self.reactants, index=range(maxnat))

        newcolumn_names = ["HOAO_a (eV)", "HOAO_b (eV)", "LUAO_a (eV)", "LUAO_b (eV)"]
        dfs = [self.HOAO_a, self.HOAO_b, self.LUAO_a, self.LUAO_b]

        for istruc in range(len(self.reactants)):
            strucid = self.reactants[istruc]
            if afo_l is None:
                mldf_original = storage.allstruc.get_value(f'{strucid}.struc_info.reactivity.afo').copy()

            else:
                mldf_original = afo_l[istruc].copy()
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

    def inter(self, spr_mode):
        identical_atom_d = find_identical_at_from_afo(self.HOAO_a,self.HOAO_b,self.LUAO_a,self.LUAO_b)
        final_react = []
        if spr_mode == 1:
            react_aa = build_react_l(self.HOAO_a, self.LUAO_a, search_mode=1)
            react_bb = build_react_l(self.HOAO_b, self.LUAO_b, search_mode=1)
            for react in [react_aa, react_bb]:
                final_react.extend(react)
        elif spr_mode == 2:
            react_ab = build_react_l(self.HOAO_a, self.LUAO_b, search_mode=1)
            react_ba = build_react_l(self.HOAO_b, self.LUAO_a, search_mode=1)
            for react in [react_ab, react_ba]:
                final_react.extend(react)

        else:
            react_aa = build_react_l(self.HOAO_a, self.LUAO_a, search_mode=1)
            react_bb = build_react_l(self.HOAO_b, self.LUAO_b, search_mode=1)
            react_ab = build_react_l(self.HOAO_a, self.LUAO_b, search_mode=1)
            react_ba = build_react_l(self.HOAO_b, self.LUAO_a, search_mode=1)

            for react in [react_aa, react_bb, react_ab, react_ba]:
                final_react.extend(react)

        final_react = sorted(refine_react_l(final_react, identical_atom_d), key=lambda x: x['gap'])
        return final_react

    def intra(self, spr_mode):
        identical_atom_d = find_identical_at_from_afo(self.HOAO_a,self.HOAO_b,self.LUAO_a,self.LUAO_b)
        final_react = []
        if spr_mode == 1:
            react_aa = build_react_l(self.HOAO_a, self.LUAO_a, search_mode=0)
            react_bb = build_react_l(self.HOAO_b, self.LUAO_b, search_mode=0)
            for react in [react_aa, react_bb]:
                final_react.extend(react)
        elif spr_mode == 2:
            react_ab = build_react_l(self.HOAO_a, self.LUAO_b, search_mode=0)
            react_ba = build_react_l(self.HOAO_b, self.LUAO_a, search_mode=0)
            for react in [react_ab, react_ba]:
                final_react.extend(react)

        else:
            react_aa = build_react_l(self.HOAO_a, self.LUAO_a, search_mode=0)
            react_bb = build_react_l(self.HOAO_b, self.LUAO_b, search_mode=0)
            react_ab = build_react_l(self.HOAO_a, self.LUAO_b, search_mode=0)
            react_ba = build_react_l(self.HOAO_b, self.LUAO_a, search_mode=0)

            for react in [react_aa, react_bb, react_ab, react_ba]:
                final_react.extend(react)

        tmp_react = sorted(refine_react_l(final_react, identical_atom_d, intra=True), key=lambda x: x['gap'])
        # del all connected bond
        allbonds = SDP(self.reactants, self.allowhydrogen, self.allowhalogen).allbonds
        final_react = []
        for item in tmp_react:
            struc = item["struc"]
            atl = item["at_combi"][0]
            trigger = True
            for item2 in allbonds:
                if item2["struc"] == struc:
                    if item2["at_combi"] == atl or len(set(atl))==1:
                        trigger=False
            if trigger:
                final_react.append(item)
        return final_react


class SDP:
    # atom index starts with 0
    def __init__(self, reactants,
                 allow_hydrogen=True,
                 allow_halogen=True,
                 bo_l=None, coord_l=None):
        self.reactants = reactants
        self.allbonds = []
        self.dfs = []
        if not bo_l is None:
            self.dfs = bo_l
            self.coord_l = coord_l
        else:
            self.dfs = [storage.allstruc.get_value(f"{_}.struc_info.reactivity.bo") for _ in self.reactants]
            self.coord_l = [storage.allstruc.get_value(f"{_}.struc_info.general.coord") for _ in self.reactants]

        for i in range(len(reactants)):
            struc = self.reactants[i]
            coord = self.coord_l[i]
            bodf = self.dfs[i]
            if bodf is not None:
                # read through df
                for index in range(bodf.shape[0]):
                    at1, at2 = bodf.iloc[index, :2]
                    at1 = int(at1) - 1
                    at2 = int(at2) - 1
                    ele1, ele2 = coord.iloc[at1, 0], coord.iloc[at2, 0]
                    if not allow_hydrogen:
                        if ele1 == "H" or ele2 == "H":
                            continue
                    if not allow_halogen:
                        if ele1 in ["F", "Cl", "cl", "Br", "I"] or ele2 in ["F", "Cl", "cl", "Br", "I"]:
                            continue
                    mbo = bodf.iloc[index, 2]
                    if mbo < 0.3:
                        continue
                    if at1 < at2:
                        at1, at2 = at2,at1
                    else:
                        at1, at2 = at1, at2
                    if bodf.shape[1] > 2:
                        embo = bodf.iloc[index, 3]
                        self.allbonds.append({"struc": struc,
                                              "at_combi": [str(at1), str(at2)],
                                              "mbo": mbo,
                                              "embo": embo})
                    else:
                        self.allbonds.append({"struc": struc,
                                              "at_combi": [str(at1), str(at2)],
                                              "mbo": mbo,
                                              })

    def build_react_l(self, mode, tol=0.01,gaplim=0.9):
        # refine react_l
        react_l = []
        for item in self.allbonds:
            struc = item["struc"]
            at1, at2 = item["at_combi"]
            mbo = item["mbo"]
            if mode == 0 and mbo>gaplim:
                continue
            if not react_l:
                react_l.append(item)
                react_l[0]["at_combi"] = [[at1, at2]]
            else:
                # loop over react_l if degenerate dissociation pair
                degenerate = False
                for item2 in react_l:
                    if int(item2["struc"]) == int(struc) and abs(mbo - item2["mbo"]) < tol:
                        degenerate = True
                        if not any([at1 in _ for _ in item2["at_combi"]]) and not any(
                                [at2 in _ for _ in item2["at_combi"]]):
                            item2["at_combi"].append([at1, at2])
                            break
                if not degenerate:
                    react_l.append(item)
                    react_l[-1]["at_combi"] = [[at1, at2]]

        if mode == 0:
            return sorted(react_l, key=lambda x: x["mbo"])
        else:
            # embo as bo
            if "embo" in react_l[0].keys():
                return sorted(react_l, key=lambda x: x["embo"])
            else:
                print("Warning: embo doesn't exist in struc info, change back to MBO assessment")
                return sorted(react_l, key=lambda x: x["mbo"])
