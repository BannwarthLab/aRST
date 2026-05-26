import sys

from aRST.script.exploration.interA import *
from aRST.script.exploration.intraA import *
from aRST.script.exploration.intraD import *
from aRST.script.head import print_simulation_details
from aRST.script.exploration.network import rxnnetwork
from aRST.script.exploration.prepare_afo import ReactivityPreparer
from aRST.script.toolbox import del_dir
from aRST.script.exploration.reactionsearch import PathEval
from itertools import combinations_with_replacement, chain, repeat


class SingleLayerSearch(PathEval):
    def __init__(self, allreactants):
        self.allreactants = allreactants
        self.prepare_reacticity()

        self.tosim_interA = self._get_all_interA_tosim_list()
        self.tosim_intraA = self._get_all_intraA_tosim_list()
        self.tosim_dissoD = self._get_all_dissoD_tosim_list()

        self.tosim_num = len(self.tosim_interA) + len(self.tosim_intraA) + len(self.tosim_dissoD)

    def prepare_reacticity(self):
        for structue in self.allreactants:
            ReactivityPreparer([structue]).ensure_AFO_ready(mute_timeit=True)
            ReactivityPreparer([structue]).ensure_weakH_ready()

    def _get_tosim_entry_list(self, react_l, num, rxntype):
        if rxntype == "dissoD":
            return [entry for entry in react_l[:num]
                    if not self.rxn_search_exist(entry)
                    ]
        else:
            return [entry for entry in react_l[:num]
                    if not self.rxn_search_exist(entry) and self.allow_H_involved_rxn(entry)
                    ]

    def _get_all_interA_tosim_list(self):
        tosim_interA = []
        ninter = storage.setting.searching_setting.inter.number
        all_struc_combi = list(combinations_with_replacement(self.allreactants, 2))
        for reactants in all_struc_combi:
            react_l = get_interA_react_l(reactants, writefile=False)
            tosim_interA_local = self._get_tosim_entry_list(react_l, ninter, "interA")
            tosim_interA.extend(tosim_interA_local)
        return tosim_interA

    def _get_all_intraA_tosim_list(self):
        tosim_intraA = []
        nintra = storage.setting.searching_setting.intra.number
        for reactant in self.allreactants:
            react_l = get_intraA_react_l([reactant], writefile=False)
            tosim_intraA_local = self._get_tosim_entry_list(react_l, nintra, "intraA")
            tosim_intraA.extend(tosim_intraA_local)
        return tosim_intraA

    def _get_all_dissoD_tosim_list(self):
        tosim_dissoD = []
        ndisso = storage.setting.searching_setting.disso.number
        for reactant in self.allreactants:
            react_l = get_intraD_react_l([reactant], writefile=False)
            react_l_filtered = [item for item in react_l if item["mbo"] < 0.9]
            tosim_dissoD_local = self._get_tosim_entry_list(react_l_filtered, ndisso, "dissoD")
            tosim_dissoD.extend(tosim_dissoD_local)
        return tosim_dissoD

    def print_progress(self, i, rc, nprint=3, failscanwd=None):
        if rc != 0:
            print(
                f"[pathsearch] {failscanwd} fails",
                flush=True,
            )
        else:
            step = max(1, self.tosim_num // nprint)
            if i % step == 0 or i == self.tosim_num:
                print(
                    f"[pathsearch] {i}/{self.tosim_num} done",
                    flush=True,
                )

    @staticmethod
    def explore_interA_entry(entry, searchcount=0, strnr="0", scanmode=0):
        if scanmode == 0:
            res, searchcount, rc = interxtbscan(entry["reactants"], entry, strnr, searchcount)

        else:
            res, searchcount, rc = intermbff(entry["reactants"], entry, strnr, searchcount)

        return res, searchcount, rc

    @staticmethod
    def explore_intraA_entry(entry, searchcount=0, strnr="0", scanmode=0):
        if scanmode == 0:
            res, searchcount, rc = intraxtbscan(entry["reactants"], entry, strnr, searchcount)

        else:
            res, searchcount, rc = intrambff(entry["reactants"], entry, strnr, searchcount)

        return res, searchcount, rc

    @staticmethod
    def explore_disso_entry(entry, searchcount=0, strnr="0", scanmode=0):
        if scanmode == 0:
            res, searchcount, rc = dissoxtbscan(entry["reactants"], entry, strnr, searchcount)

        else:
            res, searchcount, rc = dissombff(entry["reactants"], entry, strnr, searchcount)

        return res, searchcount, rc


def Cancer(reactants_ini, step=0, debugtrigger=False):
    print(f"Start doing cancer exploartions with given structues: ", reactants_ini)
    ncycle = storage.setting.searching_setting.general.cycle
    scanmode = storage.setting.searching_setting.general.scan_mode

    reactants_all = list(set(reactants_ini.copy()))
    batch_size = 50

    tseval = storage.setting.searching_setting.general.tseval
    while step < ncycle:
        print(f"begin layer {step}:")
        print()

        searchcount = 0
        print_info = []

        layersearch = SingleLayerSearch(reactants_all)
        layertasks = chain(
            zip(repeat(layersearch.explore_interA_entry), layersearch.tosim_interA),
            zip(repeat(layersearch.explore_intraA_entry), layersearch.tosim_intraA),
            zip(repeat(layersearch.explore_disso_entry), layersearch.tosim_dissoD),
        )

        # simulating all pathway
        init_rxnguess = {}
        for i, (explore_func, entry) in enumerate(layertasks, 1):
            init_rxnguess_local, searchcount, rc = explore_func(entry=entry,
                                                                searchcount=searchcount,
                                                                strnr=str(step),
                                                                scanmode=scanmode)
            init_rxnguess.update(init_rxnguess_local)


            if init_rxnguess_local:
                layersearch.print_progress(i, rc, failscanwd=init_rxnguess_local[2].get("scanwd") if rc != 0 else None)

        # eval and add new rxn info in batch
        items = list(init_rxnguess.values())
        for i in range(0, len(items), batch_size):
            batch = items[i:i + batch_size]

            for item in batch:
                reactants_info = item[0]
                endxyz_info = item[1]
                reaction_info = item[2]

                products_info, refrxnid = PathEval.followup_process_pathway(endxyz_info)
                if tseval:
                    # check if R==P
                    from collections import Counter
                    reactants  =reactants_info["reactants"]
                    reactants_items = [storage.allstruc.sid_to_strucitems[_] for _ in reactants]
                    product_items = [(" | ".join(item["mb"].split(" | ")[:storage.setting.searching_setting.general.mbcut_length]),
                                      item["charge"],
                                      item["multiplicity"]) for item in products_info]
                    bond_change = False if Counter(reactants_items) == Counter(product_items) else True
                    if not bond_change:
                        info_temp = {"reactants": reactants,
                                     "at_l": reactants_info["at_l"],
                                     "products": reactants,
                                     "scanwd": reactants_info["scanwd"],
                                     "delta_E_ts": None,
                                     "bond_change":False
                                     }
                        print_info.append(info_temp)
                        continue

                    delta_E_ts = PathEval.TS_evaluate(reactants_info, endxyz_info,debugtrigger=debugtrigger)
                    if delta_E_ts > storage.setting.searching_setting.general.elim:
                        info_temp = {"reactants": reactants_info["reactants"],
                                     "at_l": reactants_info["at_l"],
                                     "products": ["?"],
                                     "scanwd": reactants_info["scanwd"],
                                     "delta_E_ts":delta_E_ts
                                     }
                        print_info.append(info_temp)
                        continue
                info = PathEval.evalueate_path(reactants_info,
                                                 products_info,
                                                 reaction_info,
                                                 endxyz_info, refrxnid,
                                                 layer=step + 1, linkoldnode=False, delscanwd=not debugtrigger)
                print_info.append(info)
            del_dir(not debugtrigger, storage.setting.searching_setting.general.bufferwd)
            del_dir(not debugtrigger, storage.setting.searching_setting.general.scanwd)

        print_info = rxnnetwork.self_check(print_info)
        print_simulation_details(print_info, str(step))

        # screen reactants for next layer: all thermodynamic allowed generated structures
        reactants_all = list(rxnnetwork.therdy_allowed_struc)
        print("reamaining structure list: ", reactants_all, "number=", len(reactants_all))

        # # clean up buffer folder
        # clean_up_bufferwd(False, storage.setting.searching_setting.general.bufferwd)

        step += 1
