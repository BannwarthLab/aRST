import networkx as nx
import threading
import uuid
import aRST.script.initialization as storage
from collections import Counter


class ReactionNetwork:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance.graph = nx.DiGraph()
                    cls._instance._next_index = 0

                    cls._instance.rxn_idx_to_nodeid = {}
                    cls._instance.nodeid_to_rxn_idx = {}

                    cls._instance.searchpath_index = {}  # key: (sorted_Rids, atl), val=rxnid
                    cls._instance.rxn_index = {}  # key: (sorted_R_items,sorted_unreact_items, sorted_P_items), val=rxnid
        return cls._instance

    def add_newnode(self, info, edge_l=None, linkoldnode=True,registerindex=True):
        """Adds a reaction as a node and connects reactants to products.
        network[random_nodeid] =
        {info=, searchtrial=(R_sids, atl_idx), rxn=(sorted_R_items, sorted_P_items) where items=(mbcut,c,m)}

        info = {"reactants": reactants_sids,
                "unreacted": unreactants_sids,
                "at_l": atl_idx,
                "sim_uhf":sim_uhf,
                "products": products,
                "scanwd": key,
                "E_rxn": E_rxn,
                "bond_change": bond_change,
                "therdy_allow": therdy_allow,
                "layer":layer, # of which cycle step
                "reaction_exist":true or false
                "rxnid": exist if reaction_exist==0 or 1 (rxn exists)
                }"""
        nodeid = (str(uuid.uuid4().hex[:4]))
        node_idx = self._next_index
        self._next_index += 1

        self.nodeid_to_rxn_idx[nodeid] = node_idx
        self.rxn_idx_to_nodeid[node_idx] = nodeid

        reactants = info["reactants"]
        products = info["products"]
        unreact = info["unreacted"]
        at_l = info["at_l"]
        path_key = self._normal_path_pair(reactants, at_l)
        self.searchpath_index[path_key] = node_idx

        if registerindex and not "unknow" in products:
            allstruc = storage.allstruc
            R_items = [allstruc.sid_to_strucitems[i] for i in reactants]
            U_items = [allstruc.sid_to_strucitems[i] for i in unreact]
            P_items = [allstruc.sid_to_strucitems[i] for i in products]
            rxn_key = self.normal_rxn_key_from_items(R_items,U_items,P_items)
            self.rxn_index[rxn_key] = node_idx

        self.graph.add_node(nodeid, info=info)

        # link to last node
        if linkoldnode:
            if int(info["layer"]) != 0:
                lastnode = self.get_lastnode_from_node(nodeid)
                if lastnode:
                    self.graph.add_edge(lastnode, nodeid)
                if edge_l:
                    for targetnode in edge_l:
                        self.graph.add_edge(nodeid, targetnode)
        return node_idx

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
        for node in self.graph.nodes:
            info = self.get_reaction_info(node)
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
            products = info["products"]
            unreact = info["unreacted"]
            unit = (tuple(reactants), tuple(unreact), tuple(products))
            if unit in allowed_conbi:
                continue

            allowed_nodes.append(node)
            allowed_conbi.append(unit)
        return allowed_nodes

    def get_reaction_info(self, node):
        """Returns the info stored for a given reaction node."""
        if node in self.graph:
            return self.graph.nodes[node].get("info", None)
        return None

    def self_check(self, print_item_tobeupdated=None):
        import aRST.script.initialization as storage
        # update the info E_rxn since new pathway may find enenrgic lower p structure
        for node in self.graph.nodes:
            if self.graph.nodes[node].get("layer") == 0:
                continue

            info = self.get_reaction_info(node)
            products = info["products"]
            if "unknow" in products:
                continue
            if "reaction_exist" not in info.keys() or not info["reaction_exist"]:
                reactants = info["reactants"]
                E_r = [storage.allstruc.get(_).energy.sp for _ in reactants]
                E_p = [storage.allstruc.get(_).energy.sp for _ in products]
                E_rxn = sum(E_p) - sum(E_r)
                E_rxn_kcal = E_rxn * 627.503
                info["E_rxn"] = E_rxn_kcal
            else:
                # has a refid or already recalc
                refnode = self.rxn_idx_to_nodeid[info["ref_id"]]
                info["E_rxn"] = self.get_reaction_info(refnode)["E_rxn"]
        if print_item_tobeupdated:
            data = [item for item in print_item_tobeupdated if "rxnid" in item]
            for item in data:
                rxnid = item["rxnid"]
                nodeid = self.rxn_idx_to_nodeid[rxnid]
                item["E_rxn"] = self.get_reaction_info(nodeid)["E_rxn"]
            return print_item_tobeupdated

rxnnetwork = ReactionNetwork()
