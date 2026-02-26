import os
import numpy as np
import pandas as pd
from os.path import join, exists
import networkx as nx
from aRST.script.element import element
import aRST.script.initialization as storage


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
    coordinate = coord.iloc[:, 1:4].to_numpy(dtype=np.float64)
    ele_l = coord.iloc[:, 0].to_numpy()
    return coordinate, ele_l


def get_topo(coordinate, ele_l, bo):
    from aRST.script.geom import Coord
    # index starts from 1
    topo_l = []
    instance = Coord(coordinate, ele_l)
    for index in range(bo.shape[0]):  # index starts from 1
        at1, at2, mbo = bo.iloc[index, :3]
        at1 = int(at1)
        at2 = int(at2)
        dist = instance.get_ats_distance(at1 - 1, at2 - 1)
        topo_l.append({"atoms": [at1, at2], "value": dist, "mbo": float(mbo)})
    return topo_l


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
    if atcharge is not None:
        atcharge = atcharge.to_numpy().ravel()
    return atcharge


def get_strucinfo_from_path(strucid_path):
    from aRST.script.callsys import CallxTB, CallORCA, CallMolBar, CallgXTB
    from aRST.script import input_class

    coordinate, ele_l,nel,m,charge,atcharge,Esp,constrain_l,mb = (None,) * 9

    # determine coord and ele_l
    # dft opt geom > gfn2 opt geom > input geom
    if exists(join(strucid_path, "orca", "orca.xyz")):
        coordinate, ele_l = read_coord(join(strucid_path, "orca", "orca.xyz"))
    elif exists(join(strucid_path, "xtb", "xtbopt.xyz")):
        coordinate, ele_l = read_coord(join(strucid_path, "xtb", "xtbopt.xyz"))
    # find sp coord
    elif exists(join(strucid_path, "orca", "coord.xyz")):  # dft sp
        coordinate, ele_l = read_coord(join(strucid_path, "orca", "coord.xyz"))
    elif exists(join(strucid_path, "gxtb", "coord.xyz")):  # gxtb sp
        coordinate, ele_l = read_coord(join(strucid_path, "gxtb", "coord.xyz"))

    xtb = CallxTB(join(strucid_path, "xtb"))
    charge = xtb.get_charge_from_out()
    nel = int(sum([element[_] for _ in ele_l])) - charge
    atcharge = read_atcharge(join(strucid_path,"xtb", "charges"))
    Esp_gfn = xtb.get_Esp_from_out()
    # Etrv_gfn = CallxTB(join(strucid_path, "xtb", "trv")).get_Etrv()
    # Esolv_gfn = CallxTB(join(strucid_path, "xtb", "solv")).get_Esolv()
    bo = read_bo(join(strucid_path, "xtb", "wbo"))
    constrain_l = get_topo(coordinate, ele_l, bo)
    mb = CallMolBar(coordinate, ele_l).get_mb()

    if storage.setting.searching_setting.general.sp_method== "dft":
        Esp = CallORCA(join(strucid_path, "orca")).get_energy()
        m = CallORCA(join(strucid_path, "orca")).get_m()
    else:
        # can only be gxtb
        Esp=CallgXTB(join(strucid_path, "gxtb")).get_Esp_from_out()
        m = 1 if nel % 2 == 0 else 2

    if any(v is None for v in [coordinate, ele_l, nel, m, charge, atcharge, Esp, constrain_l, mb]):
        return {},False
    else:
        return {"coordinate": coordinate, "ele_l": ele_l, "nel": nel, "m": m, "c": charge, "atcharge": atcharge,
                "Esp_gfn": Esp_gfn,
                # "Etrv_gfn": Etrv_gfn,
                # "Esolv_gfn": Esolv_gfn,
                "constrain_l": constrain_l, "Esp": Esp, "mb": mb}, True

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


def write_xyz(coordinate, ele_l, name):
    pwd = os.getcwd()
    nat = len(ele_l)
    with open(f"{pwd}/{name}.xyz", "w") as f:
        f.write(f"{nat}\n\n")
        for at, (x, y, z) in zip(ele_l, coordinate):
            f.write(f"{at:2s} {x:18.14f} {y:18.14f} {z:18.14f}\n")
        f.write("\n")


class FragmentAnalyzer:
    def __init__(self, coordinate, ele_l, bo, atcharge):

        # supra struc info
        self.Snat = len(ele_l)
        self.Sbo = bo
        self.Satcharge = atcharge
        self.Scharge = np.round(np.sum(atcharge))
        self.Scoordinate = coordinate
        self.Sele_l = ele_l
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
            atcharge = self.Satcharge[list(fragment)]
            charge_flo = np.sum(atcharge)
            charge_l_flo.append(float(charge_flo))
            fragment_data.append({
                "charge_flo": charge_flo, "atcharge": atcharge})

        # Convert float fragment charges to integers while conserving total charge
        int_fragment_charges = self.adjust_fragment_charge(charge_l_flo)

        # Update fragment data with integer values
        for i, frag in enumerate(fragment_data):
            coordinate = self.Scoordinate[self.fragments[i]]
            ele_l = self.Sele_l[self.fragments[i]]
            charge_int = int(int_fragment_charges[i])
            nel = int(sum([element[_] for _ in ele_l]) - charge_int)
            frag["coordinate"] = coordinate
            frag["ele_l"] = ele_l
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

    def __init__(self, coordinate, ele_l, bo=None):

        self.nat = len(ele_l)
        self.atoms = []
        self.coordinates = []
        self.coordinates = coordinate
        self.atoms = ele_l
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
    from pydantic import BaseModel
    if exclude_keys is None:
        exclude_keys = set()

    if isinstance(d, BaseModel):
        return convert_for_json(d.model_dump(), exclude_keys)

    if isinstance(d, np.ndarray):
        return d.tolist()

    if isinstance(d, pd.DataFrame):
        return {"__dataframe__": d.to_json()}

    elif isinstance(d, list):
        return [convert_for_json(v, exclude_keys) for v in d]

    elif isinstance(d, dict):
        return {
            k: convert_for_json(v, exclude_keys)
            for k, v in d.items()
            if k not in exclude_keys
        }

    elif isinstance(d, float):
        return round(d, 8)

    return d


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
    with open(f'{storage.setting.searching_setting.general.wd}/reaction.json', "w") as file:
        file.write(json.dumps(data, indent=2, cls=CompactJSONEncoder))
        file.write('\n')


def write_strucout():
    import json
    path = f"{storage.setting.searching_setting.general.wd}/struc.json"

    with open(path, "w") as file:
        json.dump(
            convert_for_json(storage.allstruc.data, {"atomic_connection"}),
            file,
            indent=2
        )
        file.write("\n")
