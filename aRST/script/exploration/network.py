import networkx as nx
import threading
import uuid


class ReactionNetwork:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance.graph = nx.DiGraph()
                    cls._instance.nrxn_allowed = 0 # not counting layerv 0

        return cls._instance

    def add_newnode(self, info, edge_l=None,linkoldnode =True):
        """Adds a reaction as a node and connects reactants to products.
        info = {"reactants": reactants,
                "unreacted": unreactants,
                "at_l": atl,
                "sim_uhf":sim_uhf,
                "products": products,
                "scanwd": key,
                "E_rxn": E_rxn,
                "bond_changed": bondchange,
                "therdy_allow": therdy_allow,
                "layer":layer, # of which cycle step
                }"""
        nodeid = (str(uuid.uuid4()))
        self.graph.add_node(nodeid, info=info)
        if "bond_changed" in info.keys() and "therdy_allow" in info.keys():
            if info["bond_changed"] and info["therdy_allow"]:
                self.nrxn_allowed+=1

        # link to last node
        if linkoldnode:
            if int(info["layer"]) != 0:
                lastnode = self.get_lastnode_from_node(nodeid)
                if lastnode:
                    self.graph.add_edge(lastnode, nodeid)
                if edge_l:
                    for targetnode in edge_l:
                        self.graph.add_edge(nodeid, targetnode)
        return self.get_node_index_form_nodeid(nodeid)

    def add_newedge(self, node1, node2):
        self.graph.add_edge(node1, node2)

    def get_node_index_form_nodeid(self, nodeid):
        return list(self.graph.nodes).index(nodeid)

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
        if target=="0":
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
        allowed_conbi = [] # (R,U,P)
        for node in self.graph.nodes:
            info = self.get_reaction_info(node)
            # check layer
            if info["layer"] != layer:
                continue
            # check bond changed in last reaction
            if not info["bond_changed"]:
                continue
            # check thermodynamic allowed
            if not info["therdy_allow"]:
                continue
            # check if duplicated reactions
            reactants = info["reactants"]
            products = info["products"]
            unreact = info["unreacted"]
            unit =(tuple(reactants), tuple(unreact), tuple(products))
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


rxnnetwork = ReactionNetwork()
