import pickle
import os
from pathlib import Path
from functools import wraps
import aRST.script.initialization as storage
from os.path import join,exists


def load_pkl(filename):
    filename = Path(filename)
    if not filename.exists():
        return None

    with open(filename, "rb") as f:
        return pickle.load(f)

def save_pkl_atomic(data, filename):
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)

    tmp = filename.with_suffix(filename.suffix + ".tmp")
    with open(tmp, "wb") as f:
        pickle.dump(data, f)

    os.replace(tmp, filename)

def checkpoint_scan_search(dirname="checkpoint/explore"):
    check_dir = join(storage.setting.searching_setting.general.wd, dirname)
    if not exists(check_dir):
        os.makedirs(check_dir)
    def deco(func):
        @wraps(func)
        def wrapper(self, constrain_setting, rxn_type, reactantinfo, nel):
            scanwd = reactantinfo.get("scanwd")
            pkl_file = join(check_dir,f"{scanwd}.pkl")
            old_results = load_pkl(pkl_file)
            if old_results is not None:
                print(f"[restart] loaded {scanwd}")
                return old_results
            result = func(
                self,
                constrain_setting,
                rxn_type,
                reactantinfo,
                nel,
            )
            if result is not None:
                save_pkl_atomic(result, pkl_file)
                print(f"[checkpoint] saved {scanwd}")

            return result

        return wrapper

    return deco


def checkpoint_mbff_search(dirname="checkpoint/explore"):
    check_dir = join(storage.setting.searching_setting.general.wd, dirname)
    if not exists(check_dir):
        os.makedirs(check_dir)
    def deco(func):
        @wraps(func)
        def wrapper(self, atom_connection, push_l, rxn_type, reactantinfo, nel):
            scanwd = reactantinfo.get("scanwd")
            pkl_file = join(check_dir,f"{scanwd}.pkl")
            old_results = load_pkl(pkl_file)
            if old_results is not None:
                print(f"[restart] loaded {scanwd}")
                return old_results
            result = func(
                self,
                atom_connection, push_l, rxn_type, reactantinfo, nel
            )
            if result is not None:
                save_pkl_atomic(result, pkl_file)
                print(f"[checkpoint] saved {scanwd}")

            return result

        return wrapper

    return deco