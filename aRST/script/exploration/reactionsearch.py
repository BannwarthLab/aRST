from os.path import exists
from aRST.script.exploration.interA import *
from aRST.script.exploration.intraA import *
from aRST.script.exploration.intraD import *
from aRST.script.head import print_simulation_details, print_jobhead
from aRST.script.exploration.network import rxnnetwork
from aRST.script.exploration.prepare_afo import ReactivityPreparer
from aRST.script.exploration.structure import process_fragment_recursive
from aRST.script.toolbox import timeit


class SingleCycleSearch:
    def __init__(self, reactants):
        self.reactants = reactants

    @staticmethod
    def rxn_search_exist(entry):
        reactants = entry["reactants"]
        if isinstance(reactants, list):
            reactants = tuple(reactants)
        else:
            reactants = (reactants,)

        atl = entry["atl"]

        # normalize order
        if reactants != tuple(sorted(reactants, reverse=True)):
            reactants = tuple(sorted(reactants, reverse=True))
            atl = [tuple(reversed(x)) for x in atl]
        else:
            atl = [tuple(x) for x in atl]

        searchkey = (reactants, tuple(atl))
        return searchkey in rxnnetwork.searchpath_index

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
                if not self.rxn_search_exist(entry):
                    res, searchcount = interxtbscan(self.reactants, entry, strnr, searchcount)
                    output_interA.update(res)
        else:
            # call mbff
            for entry in react_l[:ninter]:
                if not self.rxn_search_exist(entry):
                    res, searchcount = intermbff(self.reactants, entry, strnr, searchcount)
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
                if not self.rxn_search_exist(entry):
                    res, searchcount = intraxtbscan(self.reactants, entry, strnr, searchcount)
                    output_intraA.update(res)
        else:
            # call mbff
            for entry in react_l[:nintra]:
                if not self.rxn_search_exist(entry):
                    res, searchcount = intrambff(self.reactants, entry, strnr, searchcount)
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
                if not self.rxn_search_exist(entry):
                    res, searchcount = dissoxtbscan(self.reactants, entry, strnr, searchcount)
                    output_D.update(res)

        else:
            # call mbff
            for entry in react_l[:ndisso]:
                if not self.rxn_search_exist(entry):
                    res, searchcount = dissombff(self.reactants, entry, strnr, searchcount)
                    output_D.update(res)
        return output_D, searchcount

    @staticmethod
    # @timeit("Separting SuparMol")
    def followup_process_pathway(endxyz_info):
        """
        output_d[f"{reactionpath}_{count}"] = [{"reactants": [struc1, struc2],
                                                    "unreacted": unreacted,
                                                    "at_l": at_l,
                                                    "simuhf": u,"scanwd":scanwd},  # reactants info


                                                   {"atomic_charge": atcharge,
                                                    "charge": charge,
                                                    "nel": totnel,
                                                    "uhf": u,
                                                    "coordinate": endxyz,
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
            "Esp": energy,
            "strucwd":rootfolder (abswd)
        }, features_p2,....]
        """
        # test if this endxyz exists
        c = endxyz_info["charge"]
        m = endxyz_info["uhf"] + 1
        mb = endxyz_info["mb"]
        refrxnid = storage.allstruc.Psupra_to_rxnid.get((mb, c, m))
        if refrxnid is not None:
            print(refrxnid, mb)
            return [],refrxnid
        else:
            bufferwd = storage.setting.searching_setting.general.bufferwd
            products_info = process_fragment_recursive(endxyz_info["coordinate"],
                                                       endxyz_info["ele_l"],
                                                       endxyz_info["charge"],
                                                       endxyz_info["nel"], bufferwd)
            return products_info,None

    @staticmethod
    def evalueate_path(reactants_info, products_info, reaction_info, endxyz_info,refrxnid=None, layer=0, linkoldnode=True):


        from collections import Counter
        reactants = reactants_info["reactants"]
        unreactants = reactants_info["unreacted"]
        simuhf = reactants_info["simuhf"]
        atl = reactants_info["at_l"]
        scanwd = reactants_info["scanwd"]
        mbcut_length = storage.setting.searching_setting.general.mbcut_length

        info_temp = {"reactants": reactants,
                     "unreacted": unreactants,
                     "at_l": atl,
                     "sim_uhf": simuhf,
                     "products": [],
                     "scanwd": scanwd,
                     # "E_rxn": E_rxn,
                     # "bond_change": bond_change,
                     # "therdy_allow": therdy_allow,
                     "layer": layer,  # of which cycle step
                     "type": reaction_info["type"]
                     # "rxnid" =
                     # reaction_exist= T/F
                     # "ref_id"
                     }

        if refrxnid is not None:
            print(info_temp["scanwd"],"!!!")
            refnodeid = rxnnetwork.rxn_idx_to_nodeid[refrxnid]
            refnodeinfo = rxnnetwork.get_reaction_info(refnodeid)
            info_temp["products"] = refnodeinfo["products"]
            info_temp["products"] = refnodeinfo["products"]
            newnode_index = rxnnetwork.add_newnode(info=info_temp, linkoldnode=linkoldnode)
            info_temp["rxnid"]=newnode_index
            return info_temp



        reactants_items = [storage.allstruc.sid_to_strucitems[_] for _ in reactants]
        unreactants_items = [storage.allstruc.sid_to_strucitems[_] for _ in unreactants]
        product_items = [(" | ".join(item["mb"].split(" | ")[:mbcut_length]),
                          item["charge"],
                          item["multiplicity"]) for item in products_info]



        search_rxn_key = rxnnetwork.normal_rxn_key_from_items(reactants_items, unreactants_items, product_items)
        refrxn_idx = rxnnetwork.rxn_index.get(search_rxn_key)
        if refrxn_idx is not None:
            refnodeid = rxnnetwork.rxn_idx_to_nodeid[refrxn_idx]
            refnodeinfo = rxnnetwork.get_reaction_info(refnodeid)
            info_temp["reaction_exist"] = True
            info_temp["ref_id"] = refrxn_idx
            info_temp["E_rxn"] = refnodeinfo["E_rxn"]  # will be updated later
            info_temp["bond_change"] = refnodeinfo["bond_change"]
            info_temp["therdy_allow"] = refnodeinfo["therdy_allow"]
            info_temp["products"] = refnodeinfo["products"]
        else:
            # new info
            info_temp["reaction_exist"] = False
            bond_change = False if Counter(reactants_items) == Counter(product_items) else True
            info_temp["bond_change"] = bond_change
            if not bond_change:
                info_temp["E_rxn"] = 0
                lastnode = rxnnetwork.get_lastnode_from_scanwd(scanwd)
                edge_l = [lastnode]
                newnode_index = rxnnetwork.add_newnode(info=info_temp, edge_l=edge_l, linkoldnode=linkoldnode)
                info_temp["rxnid"] = newnode_index
                return info_temp

            else:
                E_p = [Pfeatures["Esp"] for Pfeatures in products_info]
                E_r = [storage.allstruc.get(i).energy.sp for i in reactants]
                all_energies = E_p + E_r
                if any(e is None or not np.isfinite(e) for e in all_energies):
                    E_rxn = np.inf
                else:
                    E_rxn = sum(E_p) - sum(E_r)
                E_rxn_kcal = E_rxn * 627.503 if not np.isinf(E_rxn) else np.inf
                info_temp["E_rxn"] = E_rxn_kcal
                # print(key,E_rxn_kcal)
                if E_rxn_kcal > storage.setting.searching_setting.general.elim:
                    info_temp["therdy_allow"] = False
                    info_temp["products"] = ["unknow"]
                    newnode_index = rxnnetwork.add_newnode(info=info_temp, linkoldnode=linkoldnode)
                    info_temp["rxnid"] = newnode_index
                    return info_temp
                else:
                    info_temp["therdy_allow"] = True
        # update or add allstruc
        products = []
        for Pfeatures in products_info:
            mb = Pfeatures["mb"]
            if not "unknow" in mb:
                P_item = (
                    " | ".join(mb.split(" | ")[:mbcut_length]), Pfeatures["charge"], Pfeatures["multiplicity"])
                ref_Pid = storage.allstruc.strucitems_to_sid.get(P_item)
                if ref_Pid is not None:
                    Esp_now = Pfeatures["Esp"]
                    Esp_ref = storage.allstruc.get(ref_Pid).energy.sp
                    if Esp_now < Esp_ref:
                        storage.allstruc.update_current_structure_info(Pfeatures, ref_Pid)
                    products.append(ref_Pid)
                else:
                    sid, s = storage.allstruc.placeholde_new_sid()
                    storage.allstruc.register_new_sid(sid, s, Pfeatures, movefolder=True)
                    products.append(sid)
            else:
                products.append("unknow")
        info_temp["products"] = products
        if refrxn_idx is not None:
            newnode_index = rxnnetwork.add_newnode(info=info_temp, linkoldnode=linkoldnode, registerindex=False)
        else:
            newnode_index = rxnnetwork.add_newnode(info=info_temp, linkoldnode=linkoldnode)

            # store endxyz info
            mb = endxyz_info["mb"]
            c = endxyz_info["charge"]
            m = endxyz_info["uhf"] + 1
            storage.allstruc.register_new_SupraMol(mb, c, m, newnode_index)
        info_temp["rxnid"] = newnode_index
        return info_temp

    def workflow(self, prefix_rpid, step):
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

                products_info,refrxnid = self.followup_process_pathway(endxyz_info)
                info_l = self.evalueate_path(reactants_info,
                                             products_info,
                                             reaction_info,
                                             endxyz_info,refrxnid,
                                             layer=step + 1)
                print_info.append(info_l)

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


def ReactionSearch(reactants_ini, step=0):
    # refproducts_ini = storage.allstruc.refproducts
    print(f"Start doing reaction exploartions with given structues: ", reactants_ini)
    ncycle = storage.setting.searching_setting.general.cycle

    while step < ncycle:
        # find reaction pathes to be explored in this cycle
        if step == 0:
            print("begin cycle 0:")
            print()
            # cycle 0

            reactants_l = reactants_ini
            # add reactants_ini as 0 node
            import uuid
            nodeid = (str(uuid.uuid4().hex[:4]))
            node_idx = rxnnetwork._next_index
            rxnnetwork.graph.add_node(nodeid, nodeid=node_idx, info={"reactants": reactants_l,
                                                                     "unreacted": [],
                                                                     "products": [],
                                                                     "layer": 0})
            singlecycle = SingleCycleSearch(reactants_l)
            print_info = singlecycle.workflow(prefix_rpid="0", step=step)
            print_simulation_details(print_info, "0")
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
                print_info = singlecycle.workflow(prefix_rpid=prefix_rpid, step=step)
                print_simulation_details(print_info, nodewd)
            step += 1
    # clean up buffer_struc
    if exists(storage.setting.searching_setting.general.bufferwd):
        import shutil
        shutil.rmtree(storage.setting.searching_setting.general.bufferwd, ignore_errors=True)
