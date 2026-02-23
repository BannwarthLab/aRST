import sys
from itertools import product
import numpy as np
import pandas as pd
import aRST.script.initialization as storage


def ceil_tol(x, round_tol=0.01, tol=1e-9):
    return np.ceil((x - tol) / round_tol) * round_tol


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

def _build_atom_to_group(identical_atom_d):
    """
    identical_atom_d[struc] is a list of {"atom": [..]} groups.
    Build: atom_to_group[struc][atom] -> group_list
    """
    atom_to_group = {}
    for struc, groups in identical_atom_d.items():
        m = {}
        for g in groups:
            atoms = g["atom"]
            for a in atoms:
                m[int(a)] = [int(x) for x in atoms]
        atom_to_group[int(struc)] = m
    return atom_to_group


def refine_react_l(react_l, identical_atom_d, intra=False):
    # merge identical pais and keep the lowest gap
    # in the degenearted pairs, at can appear only once
    out = {}
    atom_to_group = _build_atom_to_group(identical_atom_d)
    for item in react_l:
        if not intra:
            s1, s2 = map(int, item["struc_combi"])
            a1, a2 = map(int, item["at_combi"])

            g1 = atom_to_group.get(s1, {}).get(a1, [a1])
            g2 = atom_to_group.get(s2, {}).get(a2, [a2])
            all_combi = list(product(g1, g2))
            if s1 == s2:
                all_combi += list(product(g2, g1))
            used_l1, used_l2 = set(), set()
            pairs = []
            for x, y in all_combi:
                if x not in used_l1 and y not in used_l2:
                    pairs.append([x, y])
                    used_l1.add(x)
                    used_l2.add(y)
            if s1 == s2:
                pairs = [[max(x, y), min(x, y)] for x, y in pairs]
            key = ("inter", s1, s2, tuple(sorted((tuple(p) for p in pairs))))
            cur_gap = item["gap"]

            if key not in out or cur_gap < out[key]["gap"]:
                out[key] = {
                    "struc_combi": [str(s1), str(s2)],
                    "at_combi": [[str(x), str(y)] for x, y in pairs],
                    "HOAO": item["HOAO"],
                    "LUAO": item["LUAO"],
                    "gap": cur_gap,
                }
        else:
            s = int(item["struc_combi"][0])
            a1, a2 = map(int, item["at_combi"])

            g1 = atom_to_group.get(s, {}).get(a1, [a1])
            g2 = atom_to_group.get(s, {}).get(a2, [a2])

            all_combi = list(product(g1, g2))

            # enforce "each atom appears once total" (your used_l list)
            used = set()
            pairs = []
            for x, y in all_combi:
                if x not in used and y not in used:
                    pairs.append([max(x, y), min(x, y)])
                    used.add(x)
                    used.add(y)

            # IMPORTANT: dedupe by each PAIR key; keep the lowest gap
            for x, y in pairs:
                key = ("intra", s, x, y)
                cur_gap = item["gap"]
                if key not in out or cur_gap < out[key]["gap"]:
                    out[key] = {
                        "struc": str(s),
                        "at_combi": [[str(x), str(y)]],
                        "HOAO": item["HOAO"],
                        "LUAO": item["LUAO"],
                        "gap": cur_gap,
                    }

    return list(out.values())

class SPR:
    # atom index starts with 0
    def __init__(self, reactants, spr_mode, afo_l=None):
        self.reactants = reactants
        # self.elculde_elements = [] if elculde_elements is None else elculde_elements

        if afo_l is None:
            maxnat = max(storage.allstruc.get_value(f'{i}.struc_info.reactivity.afo').shape[0] for i in self.reactants)
        else:
            maxnat = max(i.shape[0] for i in afo_l)

        # collecting HL values into dataframe
        self.HOAO_a = pd.DataFrame(columns=self.reactants, index=range(maxnat))
        self.HOAO_b = pd.DataFrame(columns=self.reactants, index=range(maxnat))
        self.LUAO_a = pd.DataFrame(columns=self.reactants, index=range(maxnat))
        self.LUAO_b = pd.DataFrame(columns=self.reactants, index=range(maxnat))

        # build up a elements array
        self.elements = np.full((maxnat, len(self.reactants)), np.nan)

        newcolumn_names = ["HOAO_a (eV)", "HOAO_b (eV)", "LUAO_a (eV)", "LUAO_b (eV)"]
        dfs = [self.HOAO_a, self.HOAO_b, self.LUAO_a, self.LUAO_b]

        for istruc in range(len(self.reactants)):
            strucid = self.reactants[istruc]
            if afo_l is None:
                mldf_original = storage.allstruc.get_value(f'{strucid}.struc_info.reactivity.afo').copy()
            else:
                mldf_original = afo_l[istruc].copy()

            # fill in element array
            atom_list = list(mldf_original["Atom"])
            atom_list.extend([-1] * (maxnat - len(atom_list)))
            self.elements[:, istruc] = atom_list

            for df, column_name in zip(dfs, newcolumn_names):
                tmpl = list(mldf_original.loc[:, column_name])
                for j in range(maxnat - len(tmpl)):
                    tmpl.append(np.nan)
                df.iloc[:, istruc] = tmpl  # .reset_index(drop=True)

        # Convert to numpy once
        self.HOAO_a_np = self.HOAO_a.to_numpy(dtype=float)
        self.HOAO_b_np = self.HOAO_b.to_numpy(dtype=float)
        self.LUAO_a_np = self.LUAO_a.to_numpy(dtype=float)
        self.LUAO_b_np = self.LUAO_b.to_numpy(dtype=float)

        if spr_mode == 1:
            self.pairs = [
                (self.HOAO_a_np, self.LUAO_a_np),
                (self.HOAO_b_np, self.LUAO_b_np),
            ]
        elif spr_mode == 2:
            self.pairs = [
                (self.HOAO_a_np, self.LUAO_b_np),
                (self.HOAO_b_np, self.LUAO_a_np),
            ]
        else:
            self.pairs = [
                (self.HOAO_a_np, self.LUAO_a_np),
                (self.HOAO_b_np, self.LUAO_b_np),
                (self.HOAO_a_np, self.LUAO_b_np),
                (self.HOAO_b_np, self.LUAO_a_np),
            ]

    def build_react_l(self, H, L, struc_ids, search_mode,
                      equivalent_afogap_tol=0.01,
                      exclude_elements=None,
                      exclude_pairs=None
                      ):
        """
        Output:
          list of dicts:
            {
              "struc_combi": [struc1, struc2],
              "at_combi": [at1, at2],
              "HOAO": H,
              "LUAO": L,
              "gap": gap
            }
        """
        exclude_elements = [] if exclude_elements is None else exclude_elements
        exclude_pairs = [] if exclude_pairs is None else exclude_pairs
        ex_el = {tuple(p) for p in exclude_pairs}
        ex_el |= {(b, a) for (a, b) in ex_el}
        decimals = int(round(-np.log10(equivalent_afogap_tol)))  # for print afo gap

        react_l = []
        nat, nstruc = H.shape
        for il in range(nstruc):
            Lcol = L[:, il]
            L_elements = self.elements[:, il]
            L_ok = np.isfinite(Lcol)
            if exclude_elements:
                L_ok &= ~np.isin(L_elements, exclude_elements)

            for ih in range(nstruc):
                if search_mode == 0 and il != ih:
                    continue
                if search_mode == 1 and il == ih:
                    continue
                Hcol = H[:, ih]
                H_elements = self.elements[:, ih]

                H_ok = np.isfinite(Hcol)
                if exclude_elements:
                    H_ok &= (~np.isin(H_elements, exclude_elements))

                # Broadcast to get all (iat_l, iat_h) gaps
                # gap[i_l, i_h] corresponds to L[i_l] - H[i_h]
                gap_mat = ceil_tol(Lcol[:, None] - Hcol[None, :], round_tol=equivalent_afogap_tol)
                valid = L_ok[:, None] & H_ok[None, :]
                i_l, i_h = np.where(valid)
                if i_l.size == 0:
                    continue

                if ex_el:
                    Z_l = self.elements[i_l, il]
                    Z_h = self.elements[i_h, ih]
                    keep = np.array([(int(zl), int(zh)) not in ex_el for zl, zh in zip(Z_l, Z_h)], dtype=bool)
                    if not np.any(keep):
                        continue
                    i_l = i_l[keep]
                    i_h = i_h[keep]

                gaps = gap_mat[i_l, i_h]
                hoao_vals = Hcol[i_h]
                luao_vals = Lcol[i_l]

                s_l = struc_ids[il]
                s_h = struc_ids[ih]

                if s_l <= s_h:
                    struc1, struc2 = str(s_h), str(s_l)
                    at1 = i_h
                    at2 = i_l
                else:
                    struc1, struc2 = str(s_l), str(s_h)
                    at1 = i_l
                    at2 = i_h

                for k in range(len(gaps)):
                    react_l.append({
                        "struc_combi": [struc1, struc2],
                        "at_combi": [str(int(at1[k])), str(int(at2[k]))],
                        "HOAO": float(hoao_vals[k]),
                        "LUAO": float(luao_vals[k]),
                        "gap": round(float(gaps[k]), decimals),
                    })
        return react_l

    def inter(self, equivalent_afogap_tol=0.01, exclude_elements=None, exclude_pairs=None):
        identical_atom_d = find_identical_at_from_afo(self.HOAO_a,
                                                      self.HOAO_b,
                                                      self.LUAO_a,
                                                      self.LUAO_b,
                                                      tol=equivalent_afogap_tol)
        final_react = []
        struc_ids = np.array(self.reactants, dtype=int)
        for hoao, luao in self.pairs:
            final_react.extend(
                self.build_react_l(
                    hoao,
                    luao,
                    struc_ids,
                    search_mode=1,
                    equivalent_afogap_tol=equivalent_afogap_tol,
                    exclude_elements=exclude_elements,
                    exclude_pairs=exclude_pairs
                )
            )
        final_react = sorted(refine_react_l(final_react, identical_atom_d), key=lambda x: x['gap'])
        return final_react

    def intra(self, equivalent_afogap_tol=0.01, exclude_elements=None, exclude_pairs=None):
        identical_atom_d = find_identical_at_from_afo(self.HOAO_a,
                                                      self.HOAO_b,
                                                      self.LUAO_a,
                                                      self.LUAO_b,
                                                      tol=equivalent_afogap_tol)
        final_react = []
        struc_ids = np.array(self.reactants, dtype=int)
        for hoao, luao in self.pairs:
            final_react.extend(
                self.build_react_l(
                    hoao,
                    luao,
                    struc_ids,
                    search_mode=0,
                    equivalent_afogap_tol=equivalent_afogap_tol,
                    exclude_elements=exclude_elements,
                    exclude_pairs=exclude_pairs
                )
            )

        tmp_react = sorted(refine_react_l(final_react, identical_atom_d, intra=True), key=lambda x: x['gap'])
        # del all connected bond
        allbonds = SDP(self.reactants, exclude_elements=exclude_elements, exclude_pairs=exclude_pairs).allbonds
        final_react = []
        for item in tmp_react:
            struc = item["struc"]
            atl = item["at_combi"][0]
            trigger = True
            for item2 in allbonds:
                if item2["struc"] == struc:
                    if item2["at_combi"] == atl or len(set(atl)) == 1:
                        trigger = False
            if trigger:
                final_react.append(item)
        return final_react


class SDP:
    # atom index starts with 0
    def __init__(self, reactants,
                 bo_l=None,
                 ele_ls=None,
                 exclude_elements=None, exclude_pairs=None
                 ):
        self.reactants = reactants
        self.allbonds = []
        self.dfs = []
        if bo_l is not None:
            self.dfs = bo_l
            self.ele_ls = ele_ls
        else:
            self.dfs = [storage.allstruc.get_value(f"{_}.struc_info.reactivity.bo") for _ in self.reactants]
            self.ele_ls = [storage.allstruc.get_value(f"{_}.struc_info.general.ele_l") for _ in self.reactants]
        exclude_elements = [] if exclude_elements is None else exclude_elements
        exclude_pairs = [] if exclude_pairs is None else exclude_pairs
        from aRST.script.element import element
        exclude_elements_set = set(exclude_elements)
        exclude_pairs_set = {tuple(p) for p in exclude_pairs}
        exclude_pairs_set |= {(b, a) for (a, b) in exclude_pairs_set}  # symmetric

        def keep_bond(at1, at2, mbo, z_l):
            # mbo threshold (apply to both sources)
            if mbo < 0.3:
                return None

            z1 = int(z_l[at1])
            z2 = int(z_l[at2])

            if z1 in exclude_elements_set or z2 in exclude_elements_set:
                return None
            if (z1, z2) in exclude_pairs_set:
                return None

            if at1 < at2:
                at1, at2 = at2, at1

            return at1, at2, mbo

        allbonds = []
        for i, struc in enumerate(self.reactants):
            ele_l = self.ele_ls[i]
            z_l = np.array([element[x] for x in ele_l], dtype=np.int16)

            topolist = storage.allstruc.get_value(f"{struc}.struc_info.general.atomic_connection")
            bodf = self.dfs[i]

            if topolist is not None:
                for conn in topolist:
                    at1, at2 = conn["atoms"]
                    at1 = int(at1) - 1
                    at2 = int(at2) - 1
                    mbo = float(conn["mbo"])
                    out = keep_bond(at1, at2, mbo, z_l)
                    if out is None:
                        continue
                    at1, at2, mbo = out

                    allbonds.append({"struc": struc, "at_combi": [str(at1), str(at2)], "mbo": mbo})

            else:
                if bodf is None:
                    continue

                for row in bodf.itertuples(index=False):
                    at1 = int(row[0]) - 1
                    at2 = int(row[1]) - 1
                    mbo = float(row[2])
                    out = keep_bond(at1, at2, mbo, z_l)
                    if out is None:
                        continue
                    at1, at2, mbo = out
                    allbonds.append({"struc": struc, "at_combi": [str(at1), str(at2)], "mbo": mbo})
        self.allbonds = sorted(allbonds, key=lambda x: x["mbo"])

    def build_react_l(self, tol=0.01):
        # refine react_l
        react_l = []
        for item in self.allbonds:
            struc = item["struc"]
            at1, at2 = item["at_combi"]
            mbo = item["mbo"]
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

        return sorted(react_l, key=lambda x: x["mbo"])
