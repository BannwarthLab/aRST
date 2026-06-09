import networkx as nx
import threading
import uuid
import aRST.script.initialization as storage
from collections import Counter


class ReactionNetwork:
    """"
    Each reactionpath(scanwd) is a node;


    """
    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance.graph = nx.DiGraph()
                    cls._instance._next_index = 0

                    cls._instance.scanwd_to_rxn_idx = {}
                    cls._instance.rxn_idx_to_scanwd = {} # only record first scanwd

                    cls._instance.searchpath_index = {}  # key: (sorted_Rids, atl), val=rxnid
                    cls._instance.rxn_record = {}  # key: (sorted_R_items,sorted_unreact_items, sorted_P_items), val=scanwd
                    # cls._instance.therdy_allowed_nodes = set() # record all therdy allowed nodeid
                    # cls._instance.therdy_allowed_struc = set() # record all therdy allowed generated strucres
        return cls._instance


    def add_newnode(self, info,giverxnid=True):
        """Adds a reaction as a node and connects reactants to products.

        """
        scanwd = info.get("scanwd")
        if giverxnid:
            rxnid = info.get("rxnid")
            if rxnid is None:
                rxnid = self._next_index
                info["rxnid"] = rxnid
                self._next_index += 1
            self.scanwd_to_rxn_idx[scanwd] = rxnid
            if self.rxn_idx_to_scanwd.get(rxnid) is None:
                self.rxn_idx_to_scanwd[rxnid] = scanwd

        self.graph.add_node(scanwd, info=info)

        return info



    @staticmethod
    def _normal_rxn_key(items):
        """
        items: list of (mb, charge, multiplicity)
        """
        return frozenset(Counter(items).items())

    @staticmethod
    def _normal_path_pair(R_ids, atl):
        """
           #     R_ids: sid list of R
           #     atl: acc. reactive atom pair list
           #     :return: searchtrial=(R_sids, atl_idx)
           #     """
        if len(R_ids) == 1:
            return (
                tuple(R_ids),
                tuple(tuple(x) for x in atl)
            )

        if tuple(R_ids) != tuple(sorted(R_ids, reverse=True)):
            R_ids = list(reversed(R_ids))
            atl = [tuple(reversed(x)) for x in atl]

        return (
            tuple(R_ids),
            tuple(tuple(x) for x in atl)
        )

    @classmethod
    def normal_rxn_key_from_items(cls, reactant_items,u_items, product_items):
        return (
            cls._normal_rxn_key(reactant_items),
            cls._normal_rxn_key(u_items),
            cls._normal_rxn_key(product_items),
        )



    def add_newedge(self, node1, node2):
        self.graph.add_edge(node1, node2)

    def get_layer_nodes(self, layer):
        nodes = [_ for _ in self.graph.nodes if self.get_reaction_info(_)["layer"] == int(layer)]
        return nodes

    def get_lastnode_from_node(self, node):
        info = self.get_reaction_info(node)
        if "layer" in info.keys() and info["layer"] is not None:
            layer = int(info["layer"])
        else:
            return None
        if layer == 0:
            return None
        elif layer == 1:
            return list(self.graph.nodes)[0]
        else:
            thisscanwd = self.get_reaction_info(node)["scanwd"]
            lastscanwd = "_".join(thisscanwd.split("_")[:-1])
            lastnode = [_ for _ in self.graph.nodes if "scanwd" in self.get_reaction_info(_).keys()
                        and self.get_reaction_info(_)["scanwd"] == lastscanwd][0]
            return lastnode

    def get_lastnode_from_scanwd(self, scanwd):
        if not "_" in scanwd:
            return None
        target = "_".join(scanwd.split("_")[:-1])
        if target == "0":
            return list(self.graph.nodes)[0]
        for node in self.graph.nodes:
            info = self.get_reaction_info(node)
            if "scanwd" in info.keys():
                scanwd = info["scanwd"]
                if scanwd == target:
                    return node
        return None

    def get_allowed_layer_nodes(self, layer):
        allowed_nodes = []
        allowed_conbi = []  # (R,U,P)
        badP = [None,"?","unknow","-"]
        for node in self.graph.nodes:
            info = self.get_reaction_info(node)
            products = info["products"]
            if not products:
                continue

            if any(x in badP for x in products):
                continue


            # check layer
            if info["layer"] != layer:
                continue
            # check bond changed in last reaction
            if not info["bond_change"]:
                continue
            # check thermodynamic allowed
            if not info["therdy_allow"]:
                continue
            # check if duplicated reactions
            reactants = info["reactants"]

            unreact = info["unreacted"]
            unit = (tuple(reactants), tuple(unreact), tuple(products))
            if unit in allowed_conbi:
                continue

            allowed_nodes.append(node)
            allowed_conbi.append(unit)
        return allowed_nodes

    def get_reaction_info(self, scanwd):
        """Returns the info stored for a given reaction node (scanwd)."""
        if scanwd in self.graph:
            return self.graph.nodes[scanwd].get("info", None)
        return None

    def self_check(self, print_item_tobeupdated=None):
        import aRST.script.initialization as storage
        from aRST.script.toolbox import calc_rxn_energy
        # update the info E_rxn since new pathway may find enenrgic lower p structure
        rxnid_checked = []
        for node in self.graph.nodes:
            info = self.get_reaction_info(node)
            layer = info.get("layer")

            if layer is None or layer==0:
                continue

            rxnid = info.get("rxnid")
            if rxnid is None:
                continue

            products = info["products"]
            if any(x in [None, "?", "-", "unknow"] for x in products):
                continue


            if rxnid in rxnid_checked:
                # already recalc
                ref_scanwd = self.rxn_idx_to_scanwd[rxnid]
                info["E_rxn"] = self.get_reaction_info(ref_scanwd)["E_rxn"]
                info["G_rxn"] = self.get_reaction_info(ref_scanwd)["G_rxn"]
            else:
                rxnid_checked.append(rxnid)

                reactants = info["reactants"]
                E_r = [storage.allstruc.get(_).energy.Esp_gas for _ in reactants]
                E_p = [storage.allstruc.get(_).energy.Esp_gas for _ in products]
                E_rxn_kcal = calc_rxn_energy(E_p, E_r)
                info["E_rxn"] = E_rxn_kcal

                G_r = [storage.allstruc.get(_).energy.G for _ in reactants]
                G_p = [storage.allstruc.get(_).energy.G for _ in products]
                G_rxn_kcal = calc_rxn_energy(G_p, G_r)
                info["G_rxn"] = G_rxn_kcal

        if print_item_tobeupdated:
            data = [item for item in print_item_tobeupdated if isinstance(item.get("rxnid","-"),int)]
            for item in data:
                rxnid = item["rxnid"]
                ref_scanwd= self.rxn_idx_to_scanwd[rxnid]
                item["G_rxn"] = self.get_reaction_info(ref_scanwd)["G_rxn"]
                item["E_rxn"] = self.get_reaction_info(ref_scanwd)["E_rxn"]
            return print_item_tobeupdated

rxnnetwork = ReactionNetwork()
