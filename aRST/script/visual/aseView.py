import os
from os.path import join
import json
from ase import Atoms
import pickle
import numpy as np
from ase.visualize import view
from aRST.script.toolbox import read_coord


class BaseViewer:

    @staticmethod
    def _load_json(path, error_msg):
        try:
            with open(path, "r") as f:
                return json.load(f)
        except Exception as exc:
            raise ValueError(error_msg) from exc

    @staticmethod
    def _working_dir(path):
        parent = os.path.dirname(path)
        return parent or os.getcwd()

    @staticmethod
    def _highlight_atoms(atoms, at_l):
        atoms = atoms.copy()
        tags = np.zeros(len(atoms), dtype=int)

        for pair in at_l:
            for idx in pair:
                if 0 <= idx < len(atoms):
                    tags[idx] = 1

        atoms.set_tags(tags)
        return atoms

    @classmethod
    def _make_atoms(cls, coordinates, ele_l, at_l=None):
        atoms = Atoms(
            symbols=ele_l,
            positions=np.asarray(coordinates, dtype=float),
        )

        if at_l is not None:
            atoms = cls._highlight_atoms(atoms, at_l)

        return atoms

    def _get_coordinates_and_elements(self, struc_id):
        try:
            general = self.struc_map[struc_id]["struc_info"]["general"]
            return general.get("coordinate"), general.get("ele_l")
        except KeyError as exc:
            raise ValueError(f"Structure ID {struc_id!r} not found in struc.json") from exc

    def _build_multiple_structure_atoms(self, struc_ids, offset_distance=5.0, at_l=None):
        atoms_list = []
        offset = np.zeros(3)

        for struc_id in struc_ids:
            coordinates, ele_l = self._get_coordinates_and_elements(struc_id)
            atoms = self._make_atoms(coordinates, ele_l, at_l)

            atoms.positions += offset
            atoms_list.append(atoms)

            span_x = atoms.positions[:, 0].max() - atoms.positions[:, 0].min()
            offset[0] += span_x + offset_distance

        frame = atoms_list[0]
        for atoms in atoms_list[1:]:
            frame += atoms

        return frame


class ViewSimPath(BaseViewer):
    """
    Usage:
        ViewSimPath([pathname, "struc.json", "reaction.json"]).show()
    """

    def __init__(self, inputfiles):
        if len(inputfiles) != 3:
            raise ValueError(
                "Please provide pathname, struc.json, and reaction.json. "
                "Example: aRST -v 0_0 struc.json reaction.json"
            )

        pathname, struc_path, reaction_path = inputfiles

        struc_data = self._load_json(struc_path, "Failed to load struc.json")
        reaction_data = self._load_json(reaction_path, "Failed to load reaction.json")

        self.struc_map = {item["id"]: item for item in struc_data}
        self.reaction_map = {
            item["info"]["scanwd"]: item
            for item in reaction_data
            if "info" in item and "scanwd" in item["info"]
        }

        self.pathname = pathname.strip()

        if self.pathname not in self.reaction_map:
            raise ValueError(f"{self.pathname!r} does not exist in reaction.json")

        self.wd = self._working_dir(struc_path)
        self.checkpoint_info = self._load_checkpoint()

        self._load_reaction_info()
        self.print_reaction_info()

    def _load_checkpoint(self):
        checkpoint_path = os.path.join(
            self.wd,
            "checkpoint",
            "explore",
            f"{self.pathname}.pkl",
        )

        if not os.path.exists(checkpoint_path):
            return None

        with open(checkpoint_path, "rb") as f:
            return pickle.load(f)

    def _load_reaction_info(self):
        self.rxn_info = self.reaction_map[self.pathname]["info"]

        self.reactants = self.rxn_info.get("reactants", [])
        self.products = self.rxn_info.get("products", [])
        self.at_l = self.rxn_info.get("at_l", [])
        self.rxn_type = self.rxn_info.get("type")
        self.reaction_energy = self.rxn_info.get("E_rxn")

        if len(self.reactants) > 1:
            _, first_elements = self._get_coordinates_and_elements(self.reactants[0])
            nat_1 = len(first_elements)

            self.at_l_add_nat1 = [
                [pair[0], pair[1] + nat_1]
                for pair in self.at_l
            ]
        else:
            self.at_l_add_nat1 = self.at_l

    def print_reaction_info(self):
        info = {
            "Path (scanwd)": self.pathname,
            "Type": self.rxn_type,
            "Reactants": self.reactants,
            "Simulated UHF": self.rxn_info.get("sim_uhf"),
            "Products": self.products,
            "Reactive atoms": self.at_l,
            "E_rxn": (
                f"{self.reaction_energy:.4f} kcal/mol"
                if self.reaction_energy is not None
                else None
            ),
        }

        print("=== Reaction Info ===")
        for key, value in info.items():
            print(f"{key:<15}: {value}")
        print("======================")

    def _prepare_reactants_frame(self, offset_distance=5.0):
        if len(self.reactants) == 1:
            coordinates, ele_l = self._get_coordinates_and_elements(self.reactants[0])
            return self._make_atoms(coordinates, ele_l, self.at_l)

        return self._build_multiple_structure_atoms(
            self.reactants,
            offset_distance=offset_distance,
            at_l=self.at_l_add_nat1,
        )

    def _prepare_preliminary_products_frame(self):
        xtbopt_path = os.path.join(
            self.wd,
            "scan",
            self.pathname,
            "xtbopt.xyz",
        )

        if os.path.exists(xtbopt_path):
            coordinates, ele_l = read_coord(xtbopt_path)
            return self._make_atoms(
                coordinates,
                ele_l,
                at_l=self.at_l_add_nat1,
            )

        if self.checkpoint_info is not None:
            endxyz_info = self.checkpoint_info[1]
            return self._make_atoms(
                endxyz_info["coordinate"],
                endxyz_info["ele_l"],
                at_l=self.at_l_add_nat1,
            )

        return None

    def _prepare_final_products_frame(self, offset_distance=5.0):
        return self._build_multiple_structure_atoms(
            self.products,
            offset_distance=offset_distance,
        )

    def show(self, offset_distance=5.0):
        frames = [
            self._prepare_reactants_frame(offset_distance),
            self._prepare_preliminary_products_frame(),
            self._prepare_final_products_frame(offset_distance),
        ]

        view([frame for frame in frames if frame is not None])


class ViewReactiveAtl(BaseViewer):
    def __init__(self, inputfiles):
        if len(inputfiles) != 2:
            raise ValueError(
                "Please provide react_l and struc.json. "
                "You can run `aRST --afo_features` first to obtain these files."
            )

        react_l_path, struc_path = inputfiles

        self.data = self._load_react_l(react_l_path)

        struc_data = self._load_json(struc_path, "Failed to load struc.json")
        self.struc_map = {item["id"]: item for item in struc_data}

        self.wd = self._working_dir(react_l_path)

        self.print_reaction_info()



    @staticmethod
    def _load_react_l(path):
        try:
            data = []
            with open(path, "r") as f:
                for line in f:
                    line = line.strip().rstrip(",")
                    if line:
                        data.append(json.loads(line))
            return data
        except Exception as exc:
            raise ValueError("Failed to load react_l") from exc

    def print_reaction_info(self):
        print("=== Reactive Atom Pairs ===")
        print(f"{'Rank':<6} {'atl':<15} {'mbo/gap':<15}")
        print("-" * 40)

        for rank, item in enumerate(self.data, start=1):
            atl = item.get("atl")

            if "mbo" in item:
                value = f"mbo={item['mbo']:.3f}"
            elif "gap" in item:
                value = f"gap={item['gap']:.3f} eV"
            else:
                value = "N/A"

            print(f"{rank:<6} {str(atl):<15} {value:<15}")

        print("=" * 40)

    def _prepare_atoms(self, entry, offset_distance=5.0):
        reactants = entry.get("reactants", [])
        at_l = entry.get("at_l")
        coordinates, ele_l = self._get_coordinates_and_elements(reactants[0])
        if len(reactants) == 1:
            return self._make_atoms(coordinates, ele_l, at_l), 0

        else:
            nat_1 = len(ele_l)
            return self._build_multiple_structure_atoms(
                reactants,
                offset_distance=offset_distance,
                at_l=at_l,
            ), nat_1

    def _build_frames(self, offset_distance=5.0):
        frames = []
        # get reactant atoms
        atoms =None
        nat_1=0
        for item in self.data:
            atoms, nat_1 = self._prepare_atoms(item, offset_distance=offset_distance)
            break

        for n, item in enumerate(self.data):
            atl = item.get("atl")
            at_l_add_nat1 = [
                [pair[0], pair[1] + nat_1]
                for pair in atl
            ]

            atoms = self._highlight_atoms(atoms,at_l_add_nat1)
            frames.append(atoms)
            mbo = item.get("mbo", None)
            afogap = item.get("gap", None)
            title = []

            if mbo is not None:
                title.append(f"MBO={mbo:.2f}")

            if afogap is not None:
                title.append(f"GAP={afogap:.2f} eV")

            atoms.info["name"] = " | ".join(title)

        return frames

    def show(self, offset_distance=5.0):
        frames = self._build_frames(offset_distance=offset_distance)

        view([frame for frame in frames if frame is not None])
        return