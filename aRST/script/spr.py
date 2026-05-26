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
        out[int(struc)] = []
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

            out[int(struc)].append({"atom": degenerate})
    return out


def merge_react_l(react_l, identical_atom_d, gap_tol=0.01):
    def _degenerate_group(struc, atom):
        for grp in identical_atom_d.get(struc, []):
            atoms = grp.get("atom", [])
            if atom in atoms:
                return list(atoms)
        return [atom]

    # best[key] = record with the smallest gap for that final pair key
    best = {}

    # find out beest afo val
    for item in react_l:
        s1, s2 = item["reactants"]
        a0, b0 = item["atl"]
        deg_a = _degenerate_group(s1, a0)
        deg_b = _degenerate_group(s2, b0)
        same_structure = (s1 == s2)
        if not same_structure:
            all_pairs = list(product(deg_a, deg_b))
        else:
            all_pairs = [
                (max(a, b), min(a, b))
                for a, b in product(deg_a, deg_b)
            ]

        for a, b in all_pairs:
            key = (s1, s2, a, b)
            if key not in best or item["gap"] < best[key]["gap"]:
                best[key] = {
                    "reactants": [s1, s2],
                    "atl": [[a, b]],
                    "gap": item["gap"],
                    "HOAO": item.get("HOAO"),
                    "LUAO": item.get("LUAO"),
                }

    out = list(best.values())

    def _struct_key(entry):
        s1, s2 = entry["reactants"]
        return (s1, s2, entry.get("HOAO"), entry.get("LUAO"))

    grouped = {}
    for entry in out:
        grouped.setdefault(_struct_key(entry), []).append(entry)

    merged = []
    for _, recs in grouped.items():
        recs.sort(key=lambda r: r["gap"])
        current = None
        for r in recs:
            if current is None:
                # start a new bucket
                current = {**r}
                current["atl"] = list(r["atl"])
                continue

            # merge if within tolerance of current bucket's representative gap
            if abs(r["gap"] - current["gap"]) <= gap_tol:
                current["atl"].extend(r["atl"])
            else:
                # finalize current bucket (dedupe pairs)
                seen = set()
                uniq = []
                for a, b in current["atl"]:
                    t = (a, b)
                    if t not in seen:
                        seen.add(t)
                        uniq.append([a, b])
                current["atl"] = uniq
                merged.append(current)

                current = {**r}
                current["atl"] = list(r["atl"])

        if current is not None:
            seen = set()
            uniq = []
            for a, b in current["atl"]:
                t = (a, b)
                if t not in seen:
                    seen.add(t)
                    uniq.append([a, b])
            current["atl"] = uniq
            merged.append(current)

        out = merged

    return out


class SPR:
    # atom index starts with 0
    def __init__(self, reactants, spr_mode, afo_l=None):
        self.reactants = reactants
        # self.elculde_elements = [] if elculde_elements is None else elculde_elements

        if afo_l is None:
            maxnat = max(storage.allstruc.get(i).reactivity.afo.shape[0] for i in self.reactants)
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
                mldf_original = storage.allstruc.get(strucid).reactivity.afo.copy()
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

        if np.allclose(self.HOAO_a_np, self.HOAO_b_np) and np.allclose(self.LUAO_a_np, self.LUAO_b_np):
            self.pairs = [(self.HOAO_a_np, self.LUAO_a_np)]
        elif spr_mode == 1:
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
              "reactants": [struc1, struc2],
              "atl": [at1, at2],
              "HOAO": H,
              "LUAO": L,
              "gap": gap
            }
        """
        exclude_elements = [] if exclude_elements is None else exclude_elements
        exclude_pairs = [] if exclude_pairs is None else exclude_pairs
        exclude_element_pairs = {tuple(p) for p in exclude_pairs}
        exclude_element_pairs |= {(b, a) for (a, b) in exclude_element_pairs}

        decimals = int(round(-np.log10(equivalent_afogap_tol)))

        react_l = []
        nat, nstruc = H.shape
        seen = set()
        for il in range(nstruc):
            Lcol = L[:, il]
            L_elements = self.elements[:, il]
            L_ok = np.isfinite(Lcol)
            if exclude_elements:
                L_ok &= ~np.isin(L_elements, exclude_elements)

            for ih in range(nstruc):
                if search_mode == 0 and il != ih:  # intra mode, only search gap with same col_index in HL
                    continue
                if search_mode == 1 and il == ih:  # inter mode, only search gap with diff col_index in HL (struc id can be the same!)
                    continue
                Hcol = H[:, ih]
                H_elements = self.elements[:, ih]

                H_ok = np.isfinite(Hcol)
                if exclude_elements:
                    H_ok &= (~np.isin(H_elements, exclude_elements))

                valid = L_ok[:, None] & H_ok[None, :]
                i_l, i_h = np.where(valid)
                if i_l.size == 0:
                    continue

                s_l = int(struc_ids[il])
                s_h = int(struc_ids[ih])

                # Canonicalize structure and atom ordering
                if s_l < s_h:
                    struc1, struc2 = s_h, s_l
                    at1, at2 = i_h, i_l
                    Z1 = self.elements[at1, ih]
                    Z2 = self.elements[at2, il]

                elif s_l > s_h:
                    struc1, struc2 = s_l, s_h
                    at1, at2 = i_l, i_h
                    Z1 = self.elements[at1, il]
                    Z2 = self.elements[at2, ih]

                else:
                    struc1, struc2 = s_h, s_l
                    at1 = np.maximum(i_h, i_l)
                    at2 = np.minimum(i_h, i_l)
                    Z1 = self.elements[at1, ih]
                    Z2 = self.elements[at2, ih]

                # remove excluded elements pairs
                if exclude_element_pairs:
                    keep = np.array([(int(z1), int(z2)) not in exclude_element_pairs
                                     for z1, z2 in zip(Z1, Z2)], dtype=bool)
                    if not np.any(keep):
                        continue
                    i_l = i_l[keep]
                    i_h = i_h[keep]
                    at1 = at1[keep]
                    at2 = at2[keep]

                # Broadcast to get all (iat_l, iat_h) gaps
                gaps = ceil_tol(
                    Lcol[i_l] - Hcol[i_h],
                    round_tol=equivalent_afogap_tol
                )

                hoao_vals = Hcol[i_h]
                luao_vals = Lcol[i_l]

                for k in range(len(gaps)):
                    key = (
                        struc1,
                        struc2,
                        int(at1[k]),
                        int(at2[k]),
                        float(hoao_vals[k]),
                        float(luao_vals[k]),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    react_l.append({
                        "reactants": [struc1, struc2],
                        "atl": [int(at1[k]), int(at2[k])],
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
        final_react = sorted(merge_react_l(final_react, identical_atom_d, gap_tol=equivalent_afogap_tol),
                             key=lambda x: x['gap'])
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

        tmp_react = sorted(merge_react_l(final_react, identical_atom_d, gap_tol=equivalent_afogap_tol),
                           key=lambda x: x['gap'])
        # del all connected bond & change to reactants in len(1)
        allbonds = SDP(self.reactants, exclude_elements=exclude_elements, exclude_pairs=exclude_pairs).allbonds
        final_react = []
        for item in tmp_react:
            struc = item["reactants"][0]
            atl = item["atl"][0]
            trigger = True
            for item2 in allbonds:
                if item2["reactants"] == struc:
                    if item2["atl"] == atl or len(set(atl)) == 1:
                        trigger = False
            if trigger:
                item["reactants"] = [item["reactants"][0]]
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
            self.dfs = [storage.allstruc.get(_).reactivity.bo for _ in self.reactants]
            self.ele_ls = [storage.allstruc.get(_).general.ele_l for _ in self.reactants]
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

            topolist = storage.allstruc.get(struc).general.atomic_connection
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

                    allbonds.append({"reactants": struc, "atl": [at1, at2], "mbo": mbo})

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
                    allbonds.append({"reactants": struc, "atl": [at1, at2], "mbo": mbo})
        self.allbonds = sorted(allbonds, key=lambda x: x["mbo"])

    def build_react_l(self, tol=0.01):
        react_l = []
        groups_by_struc = {}

        for item in self.allbonds:
            struc = int(item["reactants"])
            at1, at2 = item["atl"]
            mbo = item["mbo"]

            candidates = groups_by_struc.setdefault(struc, [])

            matched = False

            for group in candidates:
                if abs(mbo - group["mbo"]) < tol:
                    matched = True
                    used_atoms = group["_used_atoms"]

                    if at1 not in used_atoms and at2 not in used_atoms:
                        group["atl"].append([at1, at2])
                        used_atoms.update((at1, at2))

                    break

            if not matched:
                new_item = item.copy()
                new_item["atl"] = [[at1, at2]]
                new_item["_used_atoms"] = {at1, at2}

                candidates.append(new_item)
                react_l.append(new_item)

        for item in react_l:
            item["reactants"] = [item["reactants"]]
            item.pop("_used_atoms", None)

        return sorted(react_l, key=lambda x: x["mbo"])
