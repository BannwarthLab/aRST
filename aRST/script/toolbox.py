import os
import numpy as np
import pandas as pd
from os.path import exists
from datetime import datetime
from functools import wraps
import shutil


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

def add_noisy_to_coord(coordinate,scale=1e-4):
    # prevernt opt error results by perfert symmetry
    coordinate += np.random.normal(scale=scale, size=coordinate.shape)
    return coordinate


def get_topo(coordinate, ele_l, bo):
    topo_l = []
    for at1, at2, mbo in bo.iloc[:, :3].itertuples(index=False, name=None):
        if mbo > 0.3:
            at1 = int(at1)
            at2 = int(at2)
            dist = np.linalg.norm(coordinate[at1 - 1] - coordinate[at2 - 1])
            topo_l.append({"atoms": [at1, at2], "value": dist, "mbo": float(mbo)})
    return topo_l



def get_topo_covaradii(coordinate, ele_l,scale=1.2):
    from aRST.script.geom import r0
    from itertools import combinations
    topo_l = []
    radii_l = np.array([r0[e] for e in ele_l], dtype=float)
    pairs = list(combinations(range(1, len(ele_l) + 1), 2))
    for at1, at2 in pairs:
        dist = np.linalg.norm(coordinate[at1 - 1] - coordinate[at2 - 1])
        dist_cov = radii_l[at1 - 1] + radii_l[at2 - 1]
        bonded = dist <= scale * dist_cov
        if bonded:
            topo_l.append({"atoms": [at1, at2], "value": dist})
    return topo_l

def calc_rxn_energy(product_vals, reactant_vals):
    HARTREE_TO_KCAL = 627.503
    vals = product_vals + reactant_vals

    if any(v is None or not np.isfinite(v) for v in vals):
        return np.inf

    return (sum(product_vals) - sum(reactant_vals)) * HARTREE_TO_KCAL


def get_cm(xyffile):
    with open(xyffile, "r") as f:
        lines = f.readlines()

    if len(lines) < 2:
        raise ValueError(f"{xyffile} does not have enough lines")

    parts = lines[1].split()
    if len(parts) < 2:
        raise ValueError(f"Second line of {xyffile} must contain c and m")

    c = int(parts[0])
    m = int(parts[1])

    return c, m


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


def del_dir(trigger,wd):
    if trigger and exists(wd):
        shutil.rmtree(wd, ignore_errors=True)

def get_random_folder_name():
    import random
    import string
    random_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
    return random_name


def write_xyz(coordinate, ele_l, name, c=None, m=None):
    pwd = os.getcwd()
    nat = len(ele_l)
    cmline = "" if c is None and m is None else f"{c},{m}"
    with open(f"{pwd}/{name}.xyz", "w") as f:
        f.write(f"{nat}\n{cmline}\n")
        for at, (x, y, z) in zip(ele_l, coordinate):
            f.write(f"{at:2s} {x:18.14f} {y:18.14f} {z:18.14f}\n")
        f.write("\n")


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
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        return super().default(obj)


def write_rxnout():
    import aRST.script.initialization as storage
    from aRST.script.exploration.network import rxnnetwork
    import json

    path = f'{storage.setting.searching_setting.general.wd}/reaction.json'
    nodes = list(rxnnetwork.graph.nodes)

    with open(path, "w") as file:
        file.write("[\n")
        for i, node in enumerate(nodes):
            info = rxnnetwork.get_reaction_info(node)

            file.write("{\n")
            file.write(f'  "node": "{node}",\n')
            file.write(f'  "info": {json.dumps(info, cls=CompactJSONEncoder)}\n')
            file.write("}")

            if i != len(nodes) - 1:
                file.write(",")

            file.write("\n")
        file.write("]\n")


def to_jsonable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    else:
        return obj


def write_strucout():
    import json
    import aRST.script.initialization as storage

    path = f"{storage.setting.searching_setting.general.wd}/struc.json"
    all_sids = list(storage.allstruc.data.keys())

    with open(path, "w") as file:
        file.write("[\n")

        for i, sid in enumerate(all_sids):
            s = storage.allstruc.get(sid)
            g = s.general.model_dump(exclude={"atomic_connection"})
            e = s.energy.model_dump()

            coordinate = json.dumps(g["coordinate"], cls=CompactJSONEncoder)
            ele_l = json.dumps(g["ele_l"], cls=CompactJSONEncoder)
            atcharge = json.dumps(g["atcharge"], cls=CompactJSONEncoder)
            mb = json.dumps(g["mb"], cls=CompactJSONEncoder)
            number = json.dumps(g["number"], cls=CompactJSONEncoder)
            atomic_connection = json.dumps(
                s.general.atomic_connection,
                cls=CompactJSONEncoder
            )

            file.write("{\n")
            file.write(f'  "id": {sid},\n')
            file.write('  "struc_info": {\n')
            file.write('    "general": {\n')
            file.write(f'      "coordinate": {coordinate},\n')
            file.write(f'      "ele_l": {ele_l},\n')
            file.write(f'      "charge": {g["charge"]},\n')
            file.write(f'      "multiplicity": {g["multiplicity"]},\n')
            file.write(f'      "nel": {g["nel"]},\n')
            file.write(f'      "atcharge": {atcharge},\n')
            file.write(f'      "mb": {mb},\n')
            file.write(f'      "atomic_connection": {atomic_connection}\n')
            file.write('    },\n')
            file.write(f'    "energy": {json.dumps(e, cls=CompactJSONEncoder)}\n')
            file.write('  }\n')
            file.write('}')

            if i != len(all_sids) - 1:
                file.write(',')

            file.write('\n')

        file.write("]\n")


def timeit(label=None):
    def decorator(func):
        @wraps(func)
        def wrapper(*args,mute_timeit=False, **kwargs):
            if mute_timeit:
                return func(*args, **kwargs)
            start_time = datetime.now()

            result = func(*args, **kwargs)

            end_time = datetime.now()
            total_seconds = (end_time - start_time).total_seconds()


            print(
                f"{label} takes "
                f"{int(total_seconds // 3600)} hour "
                f"{int((total_seconds % 3600) // 60)} mins "
                f"{total_seconds % 60:.2f} sec",
                flush=True
            )

            return result

        return wrapper

    return decorator
