from os.path import exists
from aRST.script.exploration.interA import *
from aRST.script.exploration.intraA import *
from aRST.script.exploration.intraD import *
from aRST.script.head import print_simulation_details, print_jobhead
from aRST.script.exploration.network import rxnnetwork
from aRST.script.exploration.prepare_afo import ReactivityPreparer
from aRST.script.exploration.structure import process_fragment_recursive, map_fragments_topo_l_back_to_Supramol, \
    check_target_bond_change, TSstructure_NEBTS, Singlestruc
from aRST.script.toolbox import timeit, del_dir, calc_rxn_energy, get_random_folder_name
import shutil
from collections import Counter


class EntryFilter:
    @staticmethod
    def allow_H_involved_rxn(entry):
        reactants = entry["reactants"]
        atl = entry["atl"]

        structs = [storage.allstruc.get(r) for r in reactants]
        ele_ls = [s.general.ele_l for s in structs]
        weakHs = [set(s.reactivity.weakH) for s in structs]

        def is_only_single_H_structure(i):
            return len(ele_ls[i]) == 1 and ele_ls[i][0].lower() == "h"

        for a, b in atl:
            if len(reactants) == 1:
                # (struc, at)
                atoms = [(0, a), (0, b)]
            else:
                atoms = [(0, a), (1, b)]

            h_atoms = [
                (si, ai)
                for si, ai in atoms
                if ele_ls[si][ai].lower() == "h"
            ]

            # no H involved -> allow
            if not h_atoms:
                return True

            # any of structure is single-atom H ->allow
            if any(is_only_single_H_structure(si) for si, _ in h_atoms):
                return True

            # any H is strong -> not allowed
            if any(ai not in weakHs[si] for si, ai in h_atoms):
                return False

        return True

    @staticmethod
    def rxn_search_exist(entry):
        reactants = entry["reactants"]
        reactants = tuple(reactants)
        # if isinstance(reactants, list):
        #     reactants = tuple(reactants)
        # else:
        #     reactants = (reactants,)

        atl = entry["atl"]

        # normalize order
        if reactants != tuple(sorted(reactants, reverse=True)):
            reactants = tuple(sorted(reactants, reverse=True))
            atl = [tuple(reversed(x)) for x in atl]
        else:
            atl = [tuple(x) for x in atl]

        searchkey = (reactants, tuple(atl))
        return searchkey in rxnnetwork.searchpath_index


class PathEval:
    """
    Workflow for each reaction trial entry (scanwd):

    1. Get supar_endxyz from each scanend
                      │
                      │
    2.         MB item check     ─── R = P ───> terminal eval and register a rxn entry in network
                      │          ─── exist in supar_endxyz record ───> terminal eval and register entry with known info
                      │ new
    3.     (optional) TS eval    ─── fail ───> terminal eval and register a rxn entry in network
                      │
                      │ pass
       3A. barierless          3B. low-energy              3C. high-energy ───> terminal eval and register
               │                       │
               │                       │
     Pseudobarrier=5 kcal/mol          │
               │───────────────────────│
                          │
    4.       Separate endxyz and process(opt+features) each indepented fragment
                          │
                          │
    5.              MB item check ─── R = P ───> terminal eval and register a rxn entry in network
                          │
                          │       ─── at_l bond term not change ───> terminal eval and register a rxn entry in network
                          │
                          │       ─── rxn exist ───> terminal eval and register entry with known info
                          │
                          │
    6.  assign strucreID to new P and update features to exist P (may have lower energy), register rxn enrtry
                          │
                          │
    FINAL: after all rp eval finish, reaction network self check

    """

    use_g = storage.setting.searching_setting.general.evalmode == "g"
    tseval = storage.setting.searching_setting.general.tseval
    struc_energy_key = "G" if use_g else "Esp_gas"
    rxn_enenrgy_key = "G_rxn" if use_g else "E_rxn"
    ts_energy_key = "dG_ts" if use_g else "dE_ts"
    rootfolder = storage.setting.searching_setting.general.bufferwd
    elim = storage.setting.searching_setting.general.elim
    mbcut_length = storage.setting.searching_setting.general.mbcut_length
    badP = [None, "?", "-", "unknow"]

    TS_KEYS = [
        "TScoordinate",
        "TSele_l",
        "TS_Esp_gas",
        "TS_Gtrv",
        "TS_Gsolv_corr",
        "TS_G",
        "TS_barrierless",
        "dE_ts",
        "dG_ts",
    ]

    REF_COPY_KEYS = [
        "products",
        "target_change",
        "bond_change",
        "therdy_allow",
        "E_rxn",
        "G_rxn",
        "note",
        "rxnid"
    ]

    def __init__(self, reactants_info, endxyz_info, reaction_info, debugmode=False):
        self.debugmode = debugmode
        self.reactants = reactants_info["reactants"]
        self.reactants_mb = reactants_info["mb"]
        self.unreacted = reactants_info["unreacted"]
        self.at_l = reactants_info["at_l"]
        self.scanwd = reactants_info["scanwd"]
        layer = int(self.scanwd.split("_")[0]) + 1
        self.endxyz_info = endxyz_info

        self.c = endxyz_info["charge"]
        self.m = endxyz_info["uhf"] + 1
        self.coordinate = reactants_info["coordinate"]
        self.ele_l = reactants_info["ele_l"]
        self.products = []

        self.rxninfo = {
            "reactants": self.reactants,
            "unreacted": reactants_info["unreacted"],
            "at_l": self.at_l,
            "sim_uhf": endxyz_info["uhf"],

            "products": [],
            "scanwd": self.scanwd,
            "layer": layer,
            "type": reaction_info["type"],
            # "reaction_exist": None,

            "target_change": None,
            "bond_change": None,
            "therdy_allow": None,

            # TS
            "TScoordinate": None,
            "TSele_l": None,
            "TS_Esp_gas": np.inf,
            "TS_Gtrv": np.inf,
            "TS_Gsolv_corr": np.inf,
            "TS_G": np.inf,
            "TS_barrierless": None,
            "dE_ts": np.inf,
            "dG_ts": np.inf,

            "E_rxn": np.inf,
            "G_rxn": np.inf,

            "rxnid": None,
            # "ref_id": None,
            "note": "",
        }

    def rxn_exists_base_on_SupraMol(self):
        Pinit_mb = self.endxyz_info["mb"]
        Pinit_c = self.c
        Pinit_m = self.m
        ref_scanwd = storage.allstruc.Psupra_to_scanwd.get((Pinit_mb, Pinit_c, Pinit_m))
        if ref_scanwd is None:
            # register a new one
            storage.allstruc.register_new_SupraMol(Pinit_mb, Pinit_c, Pinit_m, self.scanwd)
        return ref_scanwd

    def update_info_with_ref_scanwd(self, ref_scanwd):
        current_ts_info = {
            key: self.rxninfo.get(key)
            for key in self.TS_KEYS
        }

        refnodeinfo = rxnnetwork.get_reaction_info(ref_scanwd)

        for key in self.REF_COPY_KEYS:
            if key in refnodeinfo:
                self.rxninfo[key] = refnodeinfo[key]

        # Copy reference TS by default.
        for key in self.TS_KEYS:
            if key in refnodeinfo:
                self.rxninfo[key] = refnodeinfo[key]

        # Decide which TS is better.
        current_barrier = current_ts_info.get(self.ts_energy_key, np.inf)
        ref_barrier = refnodeinfo.get(self.ts_energy_key, np.inf)

        current_has_ts = current_barrier is not None and not np.isinf(current_barrier)
        ref_has_ts = ref_barrier is not None and not np.isinf(ref_barrier)

        # If current TS is valid and better than reference TS, keep current TS.
        if current_has_ts and (not ref_has_ts or current_barrier < ref_barrier):
            for key in self.TS_KEYS:
                self.rxninfo[key] = current_ts_info[key]

        return self.rxninfo

    def check_target_bond_change(self, products_info):
        # reactants = reactants_info["reactants"]
        # atl = reactants_info["at_l"]
        map_topo_l = map_fragments_topo_l_back_to_Supramol(products_info)
        if self.rxninfo["type"] == "interA":
            struc1, struc2 = self.reactants
            nat1 = len(storage.allstruc.get(struc1).general.ele_l)
            maped_atl = [[a, b + nat1] for a, b in self.at_l]
        else:
            maped_atl = self.at_l
        target_change = check_target_bond_change(maped_atl, map_topo_l, self.rxninfo["type"])
        return target_change

    def update_rxnif_as_RequalP(self):
        self.rxninfo["note"] = "R=P"
        self.rxninfo["products"] = self.reactants
        self.rxninfo["target_change"] = False
        self.rxninfo["bond_change"] = False
        self.rxninfo["therdy_allow"] = False

        self.rxninfo["E_rxn"] = 0
        self.rxninfo["G_rxn"] = 0

        self.rxninfo["TS_Esp_gas"] = 0
        self.rxninfo["TS_Gtrv"] = 0
        self.rxninfo["TS_Gsolv_corr"] = 0
        self.rxninfo["TS_G"] = 0
        self.rxninfo["TS_barrierless"] = False
        self.rxninfo["dE_ts"] = 0
        self.rxninfo["dG_ts"] = 0

    def run(self):
        # check R = P after scan (supraMol check)
        reactants_mb = " | ".join(self.reactants_mb.split(" | ")[:self.mbcut_length])
        endxyz_mb = " | ".join(self.endxyz_info["mb"].split(" | ")[:self.mbcut_length])
        if reactants_mb == endxyz_mb:
            self.update_rxnif_as_RequalP()
            rxnnetwork.add_newnode(self.rxninfo)
            return self.rxninfo

        reactants_items = [storage.allstruc.sid_to_strucitems[_] for _ in self.reactants]
        unreactants_items = [storage.allstruc.sid_to_strucitems[_] for _ in self.unreacted]

        # check if scan-end supraMol already exist in record
        ref_scanwd = self.rxn_exists_base_on_SupraMol()
        if ref_scanwd is not None:
            # get exist P info
            products = rxnnetwork.get_reaction_info(ref_scanwd)["products"]
            if any(x in self.badP for x in products):
                rxnnetwork.add_newnode(self.rxninfo, giverxnid=False)
                return self.rxninfo
            product_items = [storage.allstruc.sid_to_strucitems[_] for _ in products]
            self.rxninfo["products"] = products

            search_rxn_key = rxnnetwork.normal_rxn_key_from_items(reactants_items, unreactants_items, product_items)
            ref_rxnidx = rxnnetwork.rxn_index.get(search_rxn_key)
            # get exist refid with known (sorted_R_items,sorted_unreact_items, sorted_P_items)
            if ref_rxnidx is not None:
                ref_scanwd = rxnnetwork.rxn_idx_to_scanwd[ref_rxnidx]
                self.update_info_with_ref_scanwd(ref_scanwd)
                self.rxninfo["note"] = "rxn exists"
                rxnnetwork.add_newnode(self.rxninfo)
                return self.rxninfo
        # defalut processing
        ## (optional) TS eval
        if self.tseval:
            TS_info = self.run_TS_eval()
            if TS_info["rc"] == 0:
                for key in TS_info:
                    if key in self.rxninfo:
                        self.rxninfo[key] = TS_info[key]
                if TS_info[self.ts_energy_key] > self.elim:
                    self.rxninfo["products"] = ["?"]
                    self.rxninfo["note"] = "high-energy TS"
                    rxnnetwork.add_newnode(self.rxninfo, giverxnid=False)
                    return self.rxninfo
            else:
                self.rxninfo["products"] = ["?"]
                self.rxninfo["note"] = "TS search fails"
                rxnnetwork.add_newnode(self.rxninfo, giverxnid=False)
                return self.rxninfo

        ## Separate endxyz and process(opt+features) each indepented fragment
        products_info = self.Splitting_supraMol(debugmode=self.debugmode)
        product_items = [(" | ".join(item["mb"].split(" | ")[:self.mbcut_length]),
                          item["charge"],
                          item["multiplicity"]) for item in products_info]

        ### check if rxn exist
        search_rxn_key = rxnnetwork.normal_rxn_key_from_items(reactants_items, unreactants_items, product_items)
        ref_rxnidx = rxnnetwork.rxn_index.get(search_rxn_key)
        if ref_rxnidx is not None:
            ref_scanwd = rxnnetwork.rxn_idx_to_scanwd[ref_rxnidx]
            self.update_info_with_ref_scanwd(ref_scanwd)
            self.rxninfo["note"] = "rxn exists"
            rxnnetwork.add_newnode(self.rxninfo)
            return self.rxninfo

        ### check if R=P
        bond_change = False if Counter(reactants_items) == Counter(product_items) else True
        self.rxninfo["bond_change"] = bond_change
        if not bond_change:
            rxnnetwork.add_newnode(self.rxninfo)
            return self.rxninfo

        ### check if desired at_l bond changed
        target_change = self.check_target_bond_change(products_info)
        self.rxninfo["target_change"] = target_change
        if not target_change:
            self.rxninfo["products"] = ['-']
            self.rxninfo["note"] = "unexpect bond-change"
            rxnnetwork.add_newnode(self.rxninfo)
            return self.rxninfo

        ### new rxn
        E_p = [Pfeatures["Esp_gas"] for Pfeatures in products_info]
        E_r = [storage.allstruc.get(i).energy.Esp_gas for i in self.reactants]
        E_rxn_kcal = calc_rxn_energy(E_p, E_r)
        self.rxninfo["E_rxn"] = E_rxn_kcal

        G_p = [Pfeatures["G"] for Pfeatures in products_info]
        G_r = [storage.allstruc.get(i).energy.G for i in self.reactants]
        G_rxn_kcal = calc_rxn_energy(G_p, G_r)
        self.rxninfo["G_rxn"] = G_rxn_kcal
        if self.use_g:
            self.rxninfo["therdy_allow"] = False if G_rxn_kcal > self.elim else True
        else:
            self.rxninfo["therdy_allow"] = False if E_rxn_kcal > self.elim else True
        self.rxninfo["note"] = "high-energy P" if not self.rxninfo["therdy_allow"] else "new rxn"

        # update or add new P
        products = []
        for Pfeatures in products_info:
            mb = Pfeatures["mb"]
            if not "unknow" in mb:
                P_item = (
                    " | ".join(mb.split(" | ")[:self.mbcut_length]), Pfeatures["charge"], Pfeatures["multiplicity"])
                ref_Pid = storage.allstruc.strucitems_to_sid.get(P_item)
                energy = Pfeatures[self.struc_energy_key]
                if ref_Pid is not None:
                    ref_energy = getattr(storage.allstruc.get(ref_Pid).energy, self.struc_energy_key)

                    if energy < ref_energy:
                        storage.allstruc.update_current_structure_info(Pfeatures, ref_Pid)
                    products.append(ref_Pid)
                else:
                    if not np.isinf(energy):
                        sid, s = storage.allstruc.placeholde_new_sid()
                        storage.allstruc.register_new_sid(sid, s, Pfeatures, movefolder=True)
                        products.append(sid)
                    else:
                        products.append("unknow")
                        self.rxninfo["note"] = "MolBar fails"
            else:
                products.append("unknow")
                self.rxninfo["note"] = "MolBar fails"
        self.rxninfo["products"] = products
        rxnnetwork.add_newnode(self.rxninfo)
        return self.rxninfo

    def run_TS_eval(self):
        info = {
            "rc": None,
            "TS_barrierless": None,

            "TScoordinate": None,
            "TSele_l": None,

            "TS_Esp_gas": np.inf,
            "TS_Gtrv": np.inf,
            "TS_Gsolv_corr": np.inf,
            "TS_G": np.inf,

            "dE_ts": np.inf,
            "dG_ts": np.inf,
        }
        # double-end TS search
        ## R-side
        Rcoordinate = self.coordinate
        Rele_l = self.ele_l

        ## P-side
        Pcoordinate_init = self.endxyz_info["coordinate"]
        Pele_l = self.endxyz_info["ele_l"]
        # relax P
        random_name = get_random_folder_name()
        comdir = join(self.rootfolder, random_name)
        create_folder(comdir)
        relaxP_features = Singlestruc(Pcoordinate_init, Pele_l, self.c, self.m).process_workflow(rootfolder=comdir,
                                                                                                 ideal=False,
                                                                                                 debugmode=self.debugmode)
        Pcoordinate = relaxP_features["coordinate"]

        # run NEBTS
        proc = TSstructure_NEBTS(Rcoordinate, Rele_l, Pcoordinate, Pele_l, self.c, self.m)
        rc = proc.run(rootfolder=self.rootfolder, debugtrigger=self.debugmode)
        info["rc"] = rc
        if rc != 0:
            return info
        TS_features = proc.features
        barrierless = TS_features.get("barrierless", False)
        info["TS_barrierless"] = barrierless
        if barrierless:
            info["dE_ts"] = 5
            info["dG_ts"] = 5
            return info

        # store TS structure if available
        info["TScoordinate"] = TS_features.get("coordinate")
        info["TSele_l"] = TS_features.get("ele_l")
        key_map = {
            "Esp_gas": "TS_Esp_gas",
            "Gtrv": "TS_Gtrv",
            "Gsolv_corr": "TS_Gsolv_corr",
            "G": "TS_G",
        }

        for src_key, dst_key in key_map.items():
            info[dst_key] = TS_features.get(src_key, np.inf)

        # calculate TS barriers
        E_TS = info["TS_Esp_gas"]
        G_TS = info["TS_G"]

        if not np.isinf(E_TS):
            E_reactants = [
                storage.allstruc.get(i).energy.Esp_gas
                for i in self.reactants
            ]
            info["dE_ts"] = calc_rxn_energy([E_TS], E_reactants)

        if not np.isinf(G_TS):
            G_reactants = [
                storage.allstruc.get(i).energy.G
                for i in self.reactants
            ]
            info["dG_ts"] = calc_rxn_energy([G_TS], G_reactants)
        return info

    def Splitting_supraMol(self, debugmode=False):
        """
                output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc1, struc2],
                                                            "unreacted": unreacted,
                                                            "at_l": at_l,
                                                            "simuhf": u,"scanwd":scanwd},  # reactants info


                                                           {"atomic_charge": atcharge,
                                                            "charge": charge,
                                                            "nel": totnel,
                                                            "uhf": u,
                                                            "c": endxyz,
                                                            "ele_l": ele_l_end,
                                                            "bo": bo,
                                                            "atomic_connection": constrain_l
                                                            "mb":mb},  # endxyz info

                                                           {"type": "interA"}  # reaction info
                                                           ]

                :return: (products_info, refrxnid)
                If the end structure already maps to an existing reaction/product entry,
                returns ([], refrxnid).

                products_info = [
                features_p1 = {
                    "coordinate": self.coordinate,
                    "ele_l": self.ele_l,
                    "nel": int(sum([element[_] for _ in self.ele_l])) - self.c,
                    "multiplicity": self.m,
                    "charge": self.c,
                    "atcharge": read_atcharge(join(rootfolder, "xtb", "charges")),
                    "mb": self.mb,
                    "atomic_connection": constrain_l,
                    "Esp_gas": Esp_gas,
                    "Gtrv":Gtrv,
                    "Gsolv_corr":Gsolv_corr,
                    "G": np.inf if any(np.isinf(x) for x in [Esp_gas, Gtrv,Gsolv_corr]) else Gtrv - Gsolv_corr + Esp_gas,
                    "strucwd":rootfolder (abswd)
                }, features_p2,....]
                """
        products_info = process_fragment_recursive(self.endxyz_info, self.rootfolder, debugmode=debugmode)
        return products_info


class SingleCycleSearch:
    def __init__(self, reactants):
        self.reactants = reactants

    @timeit("Reaction simulation")
    def explore_workflow(self, strnr="0"):
        searchcount = 0
        # output_interA, output_D, output_intraA = {}, {}, {}

        ## interA
        output_interA, searchcount = self._explore_interA(searchcount=searchcount or 0, strnr=strnr)

        ## intraA
        output_intraA, searchcount = self._explore_intraA(searchcount=searchcount or 0, strnr=strnr)

        ## dissoD
        output_D, searchcount = self._explore_disso(searchcount=searchcount or 0, strnr=strnr)

        out = {}
        for outd in [output_interA, output_D, output_intraA]:
            if outd:
                out.update(outd)
        out = [{key: value} for key, value in sorted(out.items(), key=lambda x: int(x[0].split('_')[-1]))]
        return out

    def _explore_interA(self, searchcount=0, strnr="0"):
        output_interA = {}
        intersetting = storage.setting.searching_setting.inter
        if len(self.reactants) < 2 or (intersetting.number == 0 and not intersetting.atom_pair):
            return output_interA, None
        print()
        print_jobhead("Bimolecular Association")
        from datetime import datetime
        start_time = datetime.now()
        react_l = get_interA_react_l(self.reactants)
        end_time = datetime.now()
        running_time = end_time - start_time
        total_seconds = running_time.total_seconds()
        print(f"prepare AFO gap ranking takes {int(total_seconds // 3600)} hour "
              f"{int((total_seconds % 3600) // 60)} mins {total_seconds % 60:.2f} sec", flush=True)
        scanmode = storage.setting.searching_setting.general.scan_mode
        ninter = intersetting.number
        if ninter == "all":
            ninter = len(react_l)
        print(f"Paths for Ranking 0-{ninter - 1} will be simulated: ")
        if scanmode == 0:
            # call xtb scan
            for entry in react_l[:ninter]:
                if not EntryFilter.rxn_search_exist(entry):
                    res, searchcount, rc = interxtbscan(self.reactants, entry, strnr, searchcount)
                    output_interA.update(res)
        else:
            # call mbff
            for entry in react_l[:ninter]:
                if not EntryFilter.rxn_search_exist(entry):
                    res, searchcount, rc = intermbff(self.reactants, entry, strnr, searchcount)
                    output_interA.update(res)
        return output_interA, searchcount

    def _explore_intraA(self, searchcount=0, strnr="0"):
        ## intraA
        output_intraA = {}
        intrasetting = storage.setting.searching_setting.intra
        if intrasetting.number == 0 and not intrasetting.atom_pair:
            return output_intraA, None

        print()
        print_jobhead("Intramolecular Association")
        # Testing if depotonation is energetically favorable
        if intrasetting.depotonation:
            # Testing if depotonation is energetically favorable
            from aRST.script.exploration.depotonation import updating_reacntans
            updating_reacntans(self.reactants)
            ReactivityPreparer(self.reactants).prepare_AFO(afosubwd="depotonation")
        elif intrasetting.potonation:
            # Testing if depotonation is energetically favorable
            from aRST.script.exploration.depotonation import updating_reacntans
            updating_reacntans(self.reactants, job="potonation")
            ReactivityPreparer(self.reactants).prepare_AFO(afosubwd="potonation")
        react_l = get_intraA_react_l(self.reactants)
        scanmode = storage.setting.searching_setting.general.scan_mode
        nintra = intrasetting.number
        if nintra == "all":
            nintra = len(react_l)
        print(f"Paths for Ranking 0-{nintra - 1} will be simulated: ")
        if scanmode == 0:
            # call xtb scan
            for entry in react_l[:nintra]:
                if not EntryFilter.rxn_search_exist(entry):
                    res, searchcount, rc = intraxtbscan(self.reactants, entry, strnr, searchcount)
                    output_intraA.update(res)
        else:
            # call mbff
            for entry in react_l[:nintra]:
                if not EntryFilter.rxn_search_exist(entry):
                    res, searchcount, rc = intrambff(self.reactants, entry, strnr, searchcount)
                    output_intraA.update(res)
        return output_intraA, searchcount

    def _explore_disso(self, searchcount=0, strnr="0"):
        output_D = {}
        dissosetting = storage.setting.searching_setting.disso
        if dissosetting.number == 0 and not dissosetting.atom_pair:
            return output_D, None

        print()
        print_jobhead("Dissociation")
        # Testing if depotonation is energetically favorable
        if dissosetting.depotonation:
            # Testing if depotonation is energetically favorable
            from aRST.script.exploration.depotonation import updating_reacntans
            updating_reacntans(self.reactants)
            ReactivityPreparer(self.reactants).prepare_AFO(afosubwd="depotonation")
        elif dissosetting.potonation:
            # Testing if depotonation is energetically favorable
            from aRST.script.exploration.depotonation import updating_reacntans
            updating_reacntans(self.reactants, job="potonation")
            ReactivityPreparer(self.reactants).prepare_AFO(afosubwd="potonation")
        react_l = get_intraD_react_l(self.reactants)
        scanmode = storage.setting.searching_setting.general.scan_mode
        ndisso = dissosetting.number
        if ndisso == "all":
            ndisso = len(react_l)
        print(f"Paths for Ranking 0-{ndisso - 1} will be simulated: ")
        if scanmode == 0:
            # call xtb scan
            for entry in react_l[:ndisso]:
                if not EntryFilter.rxn_search_exist(entry):
                    res, searchcount, rc = dissoxtbscan(self.reactants, entry, strnr, searchcount)
                    output_D.update(res)

        else:
            # call mbff
            for entry in react_l[:ndisso]:
                if not EntryFilter.rxn_search_exist(entry):
                    res, searchcount, rc = dissombff(self.reactants, entry, strnr, searchcount)
                    output_D.update(res)
        return output_D, searchcount

    def workflow(self, prefix_rpid, debugtrigger=False):
        from datetime import datetime
        ReactivityPreparer(self.reactants).ensure_AFO_ready()
        explore_dl = self.explore_workflow(prefix_rpid)
        print_info = []

        start_time = datetime.now()

        for output_d in explore_dl:
            for item in output_d.values():
                reactants_info = item[0]
                endxyz_info = item[1]
                reaction_info = item[2]
                rxninfo = PathEval(reactants_info, endxyz_info, reaction_info, debugtrigger).run()
                print_info.append(rxninfo)

        end_time = datetime.now()
        total_seconds = (end_time - start_time).total_seconds()

        print(
            f"EVAL RPs takes "
            f"{int(total_seconds // 3600)} hour "
            f"{int((total_seconds % 3600) // 60)} mins "
            f"{total_seconds % 60:.2f} sec",
            flush=True
        )

        # before print out, update E_r_rxn
        print_info = rxnnetwork.self_check(print_info)
        return print_info


def ReactionSearch(reactants_ini, step=0, debugtrigger=False):
    # refproducts_ini = storage.allstruc.refproducts
    print(f"Start doing reaction exploartions with given structues: ", reactants_ini)
    ncycle = storage.setting.searching_setting.general.cycle
    tseval = storage.setting.searching_setting.general.tseval
    evalmode = storage.setting.searching_setting.general.evalmode

    while step < ncycle:

        # find reaction pathes to be explored in this cycle
        if step == 0:
            print("begin cycle 0:")
            print()
            # cycle 0

            reactants_l = reactants_ini
            # # add reactants_ini as 0 node
            # import uuid
            # nodeid = (str(uuid.uuid4().hex[:4]))
            # node_idx = rxnnetwork._next_index
            # rxnnetwork.graph.add_node(nodeid, nodeid=node_idx, info={"reactants": reactants_l,
            #                                                          "unreacted": [],
            #                                                          "products": [],
            #                                                          "layer": 0})
            singlecycle = SingleCycleSearch(reactants_l)
            print_info = singlecycle.workflow(prefix_rpid="0", debugtrigger=debugtrigger)
            print_simulation_details(print_info, "0", tseval=tseval, evalmode=evalmode)
            # update reaction nexwork
            step += 1
        else:
            print("begin cycle ", step, ":")
            # cycle > 0
            nodes = rxnnetwork.get_allowed_layer_nodes(layer=step)
            nodes_scanwd = [rxnnetwork.get_reaction_info(_)["scanwd"] for _ in nodes]
            print("node = ", nodes_scanwd)
            for node, nodewd in zip(nodes, nodes_scanwd):
                print("starting exploring node ", nodewd)
                rxninfo = rxnnetwork.get_reaction_info(node)
                reactants_l = rxninfo["products"] + rxninfo["unreacted"]
                print("reactants for this node: ", reactants_l)
                prefix_rpid = rxninfo["scanwd"]
                singlecycle = SingleCycleSearch(reactants_l)
                print_info = singlecycle.workflow(prefix_rpid=prefix_rpid, debugtrigger=debugtrigger)
                print_simulation_details(print_info, nodewd, tseval=tseval, evalmode=evalmode)
            step += 1


def clean_up_bufferwd(degubtrigger, bufferwd):
    # clean up buffer_struc
    if not degubtrigger and exists(bufferwd):
        shutil.rmtree(bufferwd, ignore_errors=True)
