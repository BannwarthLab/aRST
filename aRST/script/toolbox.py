import os
import numpy as np
import pandas as pd
from os.path import join, exists
import networkx as nx
from aRST.script.element import element


class mydict(dict):
    def __init__(self, data=None):
        if data is None:
            data = {}
        super().__init__(data)

    def get_value(self, target_key, default_value=None):
        keys = target_key.split('.')
        current_data = self

        for key in keys:
            if isinstance(current_data, dict) and key in current_data:
                current_data = current_data[key]
            elif isinstance(current_data, list):
                try:
                    index = int(key)
                    current_data = current_data[index]
                except (ValueError, IndexError):
                    return current_data
            else:
                return default_value

        return current_data

    def set_value(self, target_key, new_value):
        keys = target_key.split('.')
        current_data = self

        for key in keys[:-1]:
            if isinstance(current_data, dict):
                if key not in current_data:
                    current_data[key] = {}
                current_data = current_data[key]
            else:
                raise KeyError(f"Invalid key '{key}' in the dictionary.")

        last_key = keys[-1]
        current_data[last_key] = new_value

    def compare_value(self, key1, key2):
        value1 = self.get_value(key1)
        value2 = self.get_value(key2)
        if isinstance(value1, str) and isinstance(value2, str):
            return value1 == value2
        if isinstance(value1, (int, float)) and isinstance(value2, (int, float)):
            return round(value1, 2) == round(value2, 2)
        return False


def get_all_keys(d, parent_key=""):
    keys = []

    if isinstance(d, dict):
        for k, v in d.items():
            full_key = f"{parent_key}.{k}" if parent_key else k
            keys.extend(get_all_keys(v, full_key))
    elif isinstance(d, list):
        for idx, item in enumerate(d):
            full_key = f"{parent_key}.{idx}" if parent_key else str(idx)
            keys.extend(get_all_keys(item, full_key))
    else:
        # It's a leaf value (non-dict, non-list)
        keys.append(parent_key)
    return keys


def read_coord(coordpath):
    coord = pd.read_table(coordpath, header=None, skiprows=2, sep='\s+') if exists(coordpath) else None
    return coord


def read_bo(bopath):
    bo = pd.read_csv(bopath, header=None, sep='\s+', dtype={0: int, 1: int}) \
        if exists(bopath) and os.path.getsize(bopath) != 0 else None
    return bo


def read_afo_csv(afopath):
    afo = pd.read_csv(afopath, index_col=False) if exists(afopath) else None
    if afo is not None:
        afo.columns = afo.columns.str.strip()
        afo.insert(0, 'index', range(len(afo)), allow_duplicates=False)
        afo.set_index("index")
        if 'HOAO_a (eV)' in afo.columns:
            afo = afo[['index',
                       'HOAO_a (eV)', 'HOAO_b (eV)',
                       'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
        else:
            afo = afo[['index',
                       'HOAO (eV)', 'LUAO (eV)', 'Atom']]
            afo.rename(columns={'HOAO (eV)': 'HOAO_a (eV)'})
            afo.rename(columns={'LUAO (eV)': 'LUAO_a (eV)'})
            afo.insert(len(afo.columns), "HOAO_b (eV)", afo.loc[:, "HOAO_a (eV)"].copy())
            afo.insert(len(afo.columns), "LUAO_b (eV)", afo.loc[:, "LUAO_a (eV)"].copy())
    return afo


def read_atcharge(atchargepath):
    atcharge = pd.read_table(atchargepath, header=None, sep='\s+') if exists(atchargepath) else None
    return atcharge


def get_strucinfo_from_path(strucid_path):
    from aRST.script.callsys import CallxTB, CallORCA, CallMolBar
    # determine coord
    if exists(join(strucid_path, "orca", "orca.xyz")):
        coord = read_coord(join(strucid_path, "orca", "orca.xyz"))
    elif exists(join(strucid_path, "orca", "coord.xyz")):
        coord = read_coord(join(strucid_path, "orca", "coord.xyz"))
    elif exists(join(strucid_path, "xtb", "xtbopt.xyz")):
        coord = read_coord(join(strucid_path, "orca", "xtbopt.xyz"))

    # at least it should have xtb sp res

    # elif exists(join(strucid_path, "xtb", "coord.xyz")):
    #     coord = read_coord(join(strucid_path, "orca", "coord.xyz"))
    else:
        raise ValueError(f"{strucid_path} strucinfo doesn't exist!")

    xtb = CallxTB(join(strucid_path, "xtb"))
    charge = xtb.get_charge_from_out()
    ele_l = coord.iloc[:, 0].tolist()
    nel = int(sum([element[_] for _ in ele_l])) - charge
    atcharge = read_atcharge(join(strucid_path, "charges"))
    Esp_gfn = xtb.get_Esp_from_out()
    # Etrv_gfn = CallxTB(join(strucid_path, "xtb", "trv")).get_Etrv()
    Esolv_gfn = CallxTB(join(strucid_path, "xtb", "solv")).get_Esolv()
    constrain_l, ele_l = xtb.get_topo()

    Esp_orca = CallORCA(join(strucid_path, "orca")).get_energy()
    m = CallORCA(join(strucid_path, "orca")).get_m()

    mb = CallMolBar(coord).get_mb()
    return {"coord": coord, "nel": nel, "m": m, "c": charge, "atcharge": atcharge,
            "Esp_gfn": Esp_gfn,
            # "Etrv_gfn": Etrv_gfn,
            # "Esolv_gfn": Esolv_gfn,
            "constrain_l": constrain_l, "ele_l": ele_l, "Esp_orca": Esp_orca, "mb": mb}


def create_folder(wd):
    wd = os.path.abspath(wd)
    if not exists(wd):
        os.makedirs(wd)
        os.chdir(wd)
    elif os.path.isdir(wd):
        os.chdir(wd)
    else:
        os.remove(wd)
        os.makedirs(wd)
        os.chdir(wd)
    return wd


def get_random_folder_name():
    import random
    import string
    random_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
    return random_name


def write_xyz(xyz, name):
    pwd = os.getcwd()
    # xyz.round(10)
    nat = xyz.shape[0]
    with open(f'{pwd}/{name}.xyz', 'w') as f:
        f.write(f'{nat}\n')
        f.write('\n')
        for i in range(nat):
            at, x, y, z = xyz.iloc[i, 0], xyz.iloc[i, 1], xyz.iloc[i, 2], xyz.iloc[i, 3]
            f.write('%2s %12.4f %12.4f %12.4f\n' % (at, x, y, z))
        f.write('\n')


class FragmentAnalyzer:
    def __init__(self, coord, bo, atcharge):

        # supra struc info
        self.Snat = coord.shape[0]
        self.Sbo = bo
        self.Satcharge = atcharge
        self.Scharge = np.round(atcharge.iloc[:, 0].sum())
        self.Scoord = coord
        self.Sele_l = coord.iloc[:, 0].tolist()
        self.Snel = int(sum([element[_] for _ in self.Sele_l])) - self.Scharge

        # build up graph
        self.build_cn_matrix()

        self.define_graph()
        # self.degree = [deg for node, deg in self.Graph.degree()]

        # fragment struc info
        self.fragments = [list(_) for _ in nx.connected_components(self.Graph)]

        self.fragment_data = self.get_fragment_info()

        # check if nel doesn't change
        if not self.check_nel():
            self.status = 1
        else:
            self.status = 0

    def build_cn_matrix(self):
        self.Mcn = np.zeros((self.Snat, self.Snat), dtype=int)
        # build up cn matrix base on wbo file
        if self.Sbo is not None:
            for i in range(self.Sbo.shape[0]):
                if self.Sbo.iloc[i, 2] > 0.3:
                    at1, at2 = self.Sbo.iloc[i, :2]
                    at1, at2 = int(at1), int(at2)
                    self.Mcn[at1 - 1, at2 - 1] = self.Mcn[at2 - 1, at1 - 1] = 1

    def define_graph(self):
        self.Graph = nx.Graph()
        for i in range(self.Snat):
            self.Graph.add_node(i)
            for j in range(i + 1, self.Snat):
                if self.Mcn[i, j]:
                    self.Graph.add_edge(i, j)
        G_sorted = nx.Graph()
        for node in sorted(self.Graph.nodes):
            G_sorted.add_node(node)
        for u, v, data in self.Graph.edges(data=True):
            G_sorted.add_edge(min(u, v), max(u, v), weight=data["weight"] if "weight" in data.keys() else 1)
        self.Graph = G_sorted

    def get_fragment_info(self):
        charge_l_flo = []
        fragment_data = []
        for fragment in self.fragments:
            atcharge = self.Satcharge.iloc[list(fragment)]
            charge_flo = atcharge.sum()
            charge_l_flo.append(float(charge_flo))
            fragment_data.append({
                "charge_flo": charge_flo, "atcharge": atcharge})

        # Convert float fragment charges to integers while conserving total charge
        int_fragment_charges = self.adjust_fragment_charge(charge_l_flo)

        # Update fragment data with integer values
        for i, frag in enumerate(fragment_data):
            coord = self.Scoord.iloc[self.fragments[i]]
            ele_l = coord.iloc[:, 0].tolist()
            charge_int = int(int_fragment_charges[i])
            nel = int(sum([element[_] for _ in ele_l]) - charge_int)
            frag["coord"] = coord
            frag["charge_int"] = charge_int
            frag["nel"] = nel
        return fragment_data

    def check_nel(self):
        tmpnel = sum([_["nel"] for _ in self.fragment_data])

        if tmpnel == self.Snel:
            return True
        else:
            return False

    def adjust_fragment_charge(self, charge_l_flo):

        int_charges = np.round(charge_l_flo).astype(int)  # Round to the nearest integer
        charge_difference = self.Scharge - np.sum(int_charges)  # Compute discrepancy

        if charge_difference != 0:
            residuals = np.array(charge_l_flo) - int_charges  # Find decimal parts
            adjustment_indices = np.argsort(-residuals)  # Sort by highest decimal part

            for i in range(int(abs(charge_difference))):
                idx = adjustment_indices[i]
                int_charges[idx] += np.sign(charge_difference)  # Adjust by +1 or -1
        return list(int_charges)


#
class GraphAnalyzer:

    def __init__(self, coord, bo=None):

        self.nat = coord.shape[0]
        self.atoms = []
        self.coordinates = []
        for i in range(self.nat):
            self.atoms.append(coord.iloc[i, 0])
            self.coordinates.append(np.array(coord.iloc[i, 1:4]))
        self.wbo = bo

        # get topo matrix
        self.build_cn_matrix()
        # build up graph
        self.define_edge()
        self.degree = [deg for node, deg in self.Graph.degree()]

    def build_cn_matrix(self):
        self.Mcn = np.zeros((self.nat, self.nat), dtype=int)
        if self.wbo is not None:
            # build up cn matrix base on wbo file
            for i in range(self.wbo.shape[0]):
                at1, at2 = self.wbo.iloc[i, :2]
                at1 = int(at1) - 1
                at2 = int(at2) - 1
                self.Mcn[at1, at2] = self.Mcn[at2, at1] = 1
        else:
            from aRST.script.element import r0
            # build up cn matrix base on xyz itself
            for i in range(self.nat):
                for j in range(i + 1, self.nat):
                    dist = np.linalg.norm(self.coordinates[i] - self.coordinates[j])
                    threshold = r0[self.atoms[i]] + r0[self.atoms[j]]
                    if dist <= threshold * 1.1:
                        self.Mcn[i, j] = self.Mcn[j, i] = 1

    def define_edge(self):
        self.Graph = nx.Graph()
        for i in range(self.nat):
            for j in range(i):
                if self.Mcn[i, j]:
                    bond_length = np.linalg.norm(self.coordinates[i] - self.coordinates[j])
                    self.Graph.add_edge(i, j, weight=bond_length)

        G_sorted = nx.Graph()
        for node in sorted(self.Graph.nodes):
            G_sorted.add_node(node)
        for u, v, data in self.Graph.edges(data=True):
            G_sorted.add_edge(min(u, v), max(u, v), weight=data["weight"] if "weight" in data.keys() else 1)
        self.Graph = G_sorted

    def get_t2t_path(self, tarat1: int, tarat2: int):
        from itertools import islice
        # make sure tarat1 >= tarat2
        if tarat1 < tarat2:
            tarat1, tarat2 = tarat2, tarat1
        paths = nx.shortest_simple_paths(self.Graph, source=tarat1, target=tarat2, weight='weight')
        shortest_two = list(islice(paths, 2))

        if len(shortest_two) == 1:
            d = 5  # default as 5x covalent radii
            bond_length = self.Graph[tarat1][tarat2]['weight']
            return d * bond_length
        else:
            # Straighten the ring.
            length = sum(self.Graph[u][v]['weight'] for u, v in zip(shortest_two[1][:-1], shortest_two[1][1:])) * 0.9
            return length


def convert_for_json(d, exclude_keys=None):
    """
    Recursively convert a deep nested dict to be JSON serializable.
    - Exclude specific keys.
    - Convert float values to 4 decimal places.
    - Convert DataFrame to a JSON string.
    """
    if exclude_keys is None:
        exclude_keys = set()

    if isinstance(d, pd.DataFrame):
        return {"__dataframe__": d.to_json()}  # Store as JSON string
    elif isinstance(d, list):
        return [convert_for_json(v, exclude_keys) for v in d]
    elif isinstance(d, dict):
        return {k: convert_for_json(v, exclude_keys) for k, v in d.items() if k not in exclude_keys}
    elif isinstance(d, float):
        return round(d, 8)  # Round float values to 8 decimal places
    return d  # int, str, etc., remain unchanged


import json


class CompactJSONEncoder(json.JSONEncoder):
    def iterencode(self, obj, _one_shot=False):
        """Override default encoding to keep lists compact."""
        if isinstance(obj, list):
            return "[" + ",".join(map(self.encode, obj)) + "]"
        return super().iterencode(obj, _one_shot)


def write_rxnout():
    import aRST.script.initialization as storage
    from aRST.script.exploration.network import rxnnetwork
    import json
    data = {
        "nodes": [
            {
                "node": node,
                "info": rxnnetwork.get_reaction_info(node)}
            for node in rxnnetwork.graph.nodes
        ],
        "edges": list(rxnnetwork.graph.edges)
    }
    with open(f'{storage.setting.wd0}/reaction.json', "w") as file:
        file.write(json.dumps(data, indent=2, cls=CompactJSONEncoder))
        file.write('\n')


def write_strucout():
    import aRST.script.initialization as storage
    import json
    with open(f'{storage.setting.wd0}/struc.json', "w") as file:
        exclude_keys = ["atomic_connection"]
        file.write(json.dumps(convert_for_json(storage.allstruc.data, exclude_keys), indent=2))
        file.write('\n')
