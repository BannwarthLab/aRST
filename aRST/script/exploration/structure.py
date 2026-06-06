import os
import sys
import warnings
from os.path import exists, join
import numpy as np
from aRST.script.callsys import CallxTB, CallORCA, CallgXTB, CallMolBar, CallFFiTS
import aRST.script.initialization as storage
from aRST.script.toolbox import create_folder, read_coord, read_atcharge, write_xyz, read_bo, get_topo, del_dir, \
    get_random_folder_name, add_noisy_to_coord
import threading
from aRST.script import structure_class
import networkx as nx
from aRST.script.element import element
import shutil


class Singlestruc:
    # generate electronic properties of a given single structure
    def __init__(self, coordinate, ele_l, c, m, mb=None):
        self.coordinate = coordinate
        self.ele_l = ele_l
        self.c = c
        self.m = m
        if mb is None:
            self.mb = CallMolBar(coordinate, ele_l).get_mb()
        else:
            self.mb = mb

    def structure_exist(self):
        strucitems_to_sid = storage.allstruc.strucitems_to_sid
        mb_cut = " | ".join(self.mb.split(" | ")[:storage.setting.searching_setting.general.mbcut_length])
        sid = strucitems_to_sid.get((mb_cut, self.c, self.m))
        if sid:
            return True, sid
        else:
            return False, None

    def _prepare_workdir_and_write_coord(self, comdir):
        if not exists(comdir):
            create_folder(comdir)
        else:
            os.chdir(comdir)
        write_xyz(self.coordinate, self.ele_l, "coord", c=self.c, m=self.m)

    def _run_opt(self, coordinate_in, comdir, engine="xtb"):
        create_folder(comdir)
        write_xyz(coordinate_in, self.ele_l, "coord", c=self.c, m=self.m)
        if engine == "orca":
            rc = CallORCA(comdir, custom_settings={"charge": self.c, "m": self.m, "caltyp": "SP OPT"}).call(
                type="opt")
            if rc == 0:
                coordinate_out, self.ele_l = read_coord(join(comdir, "orca.xyz"))
            else:
                coordinate_out = coordinate_in
                print(f'Warning: struc in {comdir} opt at ORCA failed, continue with unopt struc')
        elif engine == "xtb":
            rc = CallxTB(comdir, custom_settings={"charge": self.c, "uhf": self.m - 1}).opt()
            if rc == 0:
                coordinate_out, self.ele_l = read_coord(join(comdir, "xtbopt.xyz"))
            else:
                coordinate_out = coordinate_in
                print(f'Warning: struc in {comdir} opt at xtb failed, continue with unopt struc')
        else:
            raise ValueError(f"Unknown engine for opt: {engine}")
        # update the mb
        if rc==0:
            mb = CallMolBar(coordinate_out, self.ele_l).get_mb()
            if self.mb != mb:
                # print(f'Warning: Molbar code of struc in {comdir} changed after opt')
                self.mb = mb
        return coordinate_out, self.ele_l, rc

    def _run_trv(self, comdir):
        self._prepare_workdir_and_write_coord(comdir)
        proc = CallORCA(comdir, custom_settings={
            "alpb": storage.setting.xtb_setting.alpb,
            "charge": self.c,
            "m": self.m,
        })
        Gtrv = np.inf
        rc = proc.call_trv()
        if rc == 0:
            Gtrv = proc.get_Gtrv_from_out()
        return Gtrv, rc

    def _run_solv_corr(self, comdir):
        self._prepare_workdir_and_write_coord(comdir)
        proc = CallxTB(comdir, custom_settings={"charge": self.c,
                                                "uhf": self.m - 1,
                                                })
        rc_gas, rc_solv = proc.solv_correct()
        Gsolv_corr = np.inf
        if rc_solv == 0 and rc_gas == 0:
            Gsolv_corr = proc.get_Gsolv_corr_from_out()
            rc = 0
        else:
            rc = 1
        return Gsolv_corr, rc

    def _run_sp(self, comdir, engine="gxtb", phase="gas"):
        # gas E_sp
        self._prepare_workdir_and_write_coord(comdir)
        if engine == "gxtb":
            proc = CallgXTB(comdir, custom_settings={"charge": self.c,
                                                     "uhf": self.m - 1})
            # write .charg and .uhf file
            with open(".CHRG", "w") as file:
                file.write(f"{self.c}\n")
            with open(".UHF", "w") as file:
                file.write(f"{self.m - 1}\n")
            rc = proc.sp()
            if rc == 0:
                energy = proc.get_Esp_from_out()
            else:
                energy = np.inf
                print(f'Warning: struc in {comdir} sp at gxtb failed, continue anyhow')

        elif engine == "orca":
            proc = CallORCA(comdir, custom_settings={"charge": self.c,
                                                     "m": self.m,
                                                     "caltyp": "SP",
                                                     "cpcm": None if phase == "gas" else storage.setting.orca_setting.cpcm})
            rc = proc.call(type="sp")
            if rc == 0:
                energy = proc.get_energy_from_out()
            else:
                energy = np.inf
                print(f'Warning: struc in {comdir} sp at ORCA failed, continue anyhow')
        else:
            raise ValueError(f"Unknown engine for sp: {engine}")
        if energy is None:
            energy = np.inf
            print(f'Warning: struc in {comdir} sp prop fails, structure energy will be set as inf')
        return energy, rc

    def _idealize(self, topo_l):
        from molbar.idealize import idealize_structure_from_coordinates
        # Define constraints
        extra_constraints = {
            'bond_order_assignment': True,
            'cycle_detection': True,
            'set_edges': False,
            'set_angles': True,
            'set_dihedrals': False,
            'set_repulsion': True,
            'repulsion_charge': 100,
            'constraints': {'bonds': topo_l}
        }

        # Call the idealization function
        res = idealize_structure_from_coordinates(self.coordinate, self.ele_l,
                                                  input_constraint=extra_constraints)
        coordinate, ele_l = res[1], res[2]
        if coordinate is not None and ele_l is not None:
            self.coordinate = coordinate
            self.ele_l = ele_l
        return coordinate, ele_l

    def process_workflow(self, rootfolder, ideal=True, opt=True, debugmode=False):

        # ideal->opt/sp -> collect res

        from aRST.script.toolbox import get_topo_covaradii
        if len(self.ele_l) > 1 and ideal:
            topo_l_simple = get_topo_covaradii(self.coordinate, self.ele_l)
            G = StrucGraph(self.coordinate, self.ele_l, topo_l_simple).graph
            if nx.is_connected(G):
                self._idealize(topo_l_simple)

        ## opt (defalut with solvation model)
        rc_opt = 0
        if opt:
            # add a tiny noisy to coordinante first
            coordinate_in = add_noisy_to_coord(self.coordinate,scale=1e-3)
            dftopt = storage.setting.read.general.dftopt or storage.setting.orca_setting.dftopt
            if len(self.ele_l) > 1:
                # do gfn2 opt
                comdir = join(rootfolder, "xtb")
                # coordinate_in, coordinate_out, comdir, engine="xtb"
                res = self._run_opt(coordinate_in,comdir, engine="xtb")
                rc_opt = res[2]
                if rc_opt==0:
                    coordinate_out = res[0]
                    if dftopt:
                        # further DFT opt
                        coordinate_in = coordinate_out
                        comdir = join(rootfolder, "orca_opt")
                        res = self._run_opt(coordinate_in, comdir, engine="orca")
                        rc_opt = res[2]
                        coordinate_out = res[0] if rc_opt==0 else coordinate_in
                    if rc_opt==0:
                        self.coordinate = coordinate_out

        ## sp at gas phase
        if storage.setting.searching_setting.general.sp_method == "dft":
            comdir = join(rootfolder, "orca_sp")
            Esp_gas, rc_sp = self._run_sp(comdir, engine="orca", phase="gas")
        elif storage.setting.searching_setting.general.sp_method == "gxtb":
            comdir = join(rootfolder, "gxtb")
            Esp_gas, rc_sp = self._run_sp(comdir, engine="gxtb", phase="gas")
        else:
            raise ValueError("sp_method keyword wrong")

        ## hess and solv calc when eval mode is free gibbs energy
        Gtrv = np.inf
        Gsolv_corr = np.inf
        rc_trv, rc_solv = 0, 0
        if storage.setting.searching_setting.general.evalmode == "g":
            ### hess calc
            comdir = join(rootfolder, "orca_hess")
            Gtrv, rc_trv = self._run_trv(comdir)

            ### solv
            comdir = join(rootfolder, "solv_corr")
            Gsolv_corr, rc_solv = self._run_solv_corr(comdir)  # E_solv - E_gas

        # collect data
        if len(self.ele_l) > 1:
            bo = read_bo(join(rootfolder, "xtb", "wbo"))
            constrain_l = get_topo(self.coordinate, self.ele_l, bo)
            atcharge = read_atcharge(join(rootfolder, "xtb", "charges"))
        else:
            constrain_l = []
            atcharge = [self.c]
        features = {
            "coordinate": self.coordinate,
            "ele_l": self.ele_l,
            "nel": int(sum([element[_] for _ in self.ele_l])) - self.c,
            "multiplicity": self.m,
            "charge": self.c,
            "atcharge": atcharge,
            "mb": self.mb,
            "atomic_connection": constrain_l,
            "Esp_gas": Esp_gas,
            "Gtrv": Gtrv,
            "Gsolv_corr": Gsolv_corr,
            "G": np.inf if any(np.isinf(x) for x in [Esp_gas, Gtrv, Gsolv_corr]) else Gtrv - Gsolv_corr + Esp_gas,
            "strucwd": rootfolder
        }
        # if debugmode and any(x == 1 for x in [rc_opt, rc_sp, rc_trv, rc_solv]):
        # else:
        #     del_dir(not debugmode, rootfolder)

        return features

    @staticmethod
    def read_features_from_struc_folder(rootfolder, charge, m):
        # update coord and ele_l: dftopt>gxtb>xtb

        if exists(join(rootfolder, "orca_opt", "orca.xyz")):
            coordinate, ele_l = read_coord(join(rootfolder, "orca_opt", "orca.xyz"))


        elif exists(join(rootfolder, "xtb", "xtbopt.xyz")):
            coordinate, ele_l = read_coord(join(rootfolder, "xtb", "xtbopt.xyz"))


        elif exists(join(rootfolder, "orca_sp", "coord.xyz")):  # dft sp
            coordinate, ele_l = read_coord(join(rootfolder, "orca_opt", "coord.xyz"))


        elif exists(join(rootfolder, "gxtb", "coord.xyz")):  # gxtb sp
            coordinate, ele_l = read_coord(join(rootfolder, "gxtb", "coord.xyz"))


        else:
            return {}, False

        if storage.setting.searching_setting.general.sp_method == "dft":
            Esp_gas = CallORCA(join(rootfolder, "orca_sp")).get_energy_from_out()
            if m is None:
                m = CallORCA(join(rootfolder, "orca_sp")).get_m()
        else:
            Esp_gas = CallgXTB(join(rootfolder, "gxtb")).get_Esp_from_out()

        Gtrv = np.inf
        Gsolv_corr = np.inf
        if storage.setting.searching_setting.general.evalmode == "g":
            Gtrv = CallORCA(join(rootfolder, "orca_hess")).get_Gtrv_from_out()
            Gsolv_corr = CallxTB(join(rootfolder, "solv_corr")).get_Gsolv_corr_from_out()

        bo = read_bo(join(rootfolder, "xtb", "wbo"))
        constrain_l = get_topo(coordinate, ele_l, bo)
        atcharge = read_atcharge(join(rootfolder, "xtb", "charges"))
        mb = CallMolBar(coordinate, ele_l).get_mb()
        features = {
            "coordinate": coordinate,
            "ele_l": ele_l,
            "nel": int(sum([element[_] for _ in ele_l])) - charge,
            "multiplicity": m,
            "charge": charge,
            "atcharge": atcharge,
            "mb": mb,
            "atomic_connection": constrain_l,
            "Esp_gas": Esp_gas,
            "Gtrv": Gtrv,
            "Gsolv_corr": Gsolv_corr,
            "G": np.inf if any(np.isinf(x) for x in [Esp_gas, Gtrv, Gsolv_corr]) else Gtrv - Gsolv_corr + Esp_gas,
            "strucwd": rootfolder
        }
        if any(v is None for v in features.values()):
            return {}, False

        return features, True


class Supramolecule:
    # analyse the possible Supramolecule
    def __init__(self, Supra_coordinate, Supra_ele_l, Supra_atcharge, Supra_topo_l, totc, totnel,
                 global_atomidx=None, ):
        self.Supra_coordinate = Supra_coordinate
        self.Supra_ele_l = Supra_ele_l
        self.Supra_atcharge = Supra_atcharge
        self.Supra_topo_l = Supra_topo_l
        self.totc = totc
        self.totnel = totnel
        if global_atomidx is None:
            global_atomidx = list(range(1, len(Supra_ele_l) + 1))
        self.global_atomidx = global_atomidx
        self.Sgraph = StrucGraph(coordinate=Supra_coordinate, ele_l=Supra_ele_l, topo_l=Supra_topo_l).graph

        self.Fragments_coordinate_list = []
        self.Fragments_ele_l_list = []
        self.Fragments_c_list = []
        self.Fragments_atcharge_list = []
        self.Fragments_nel_list = []
        self.Fragments_topo_list = []
        self.Fragments_global_atomidx_list = []

    def split(self):
        components = list(nx.connected_components(self.Sgraph))

        frag_coord_list = []
        frag_ele_list = []
        frag_atcharge_list = []
        frag_c_float_list = []
        frag_nel_list = []
        frag_topo_list = []
        frag_global_atomidx_list = []

        # print("sum Supra_atcharge =", float(np.sum(self.Supra_atcharge)))
        # print("totc =", self.totc)
        for comp in components:

            comp_sorted = sorted(comp)
            atom_idxs0 = [i - 1 for i in comp_sorted]
            frag_global_atomidx = [
                self.global_atomidx[i - 1]
                for i in comp_sorted
            ]

            # map supra global 1-based -> fragment local 1-based
            global_to_local = {g: l for l, g in enumerate(comp_sorted, start=1)}

            frag_coord = self.Supra_coordinate[atom_idxs0]
            frag_ele = [self.Supra_ele_l[i] for i in atom_idxs0]
            frag_atcharge = [self.Supra_atcharge[i] for i in atom_idxs0]

            frag_c_float = float(np.sum(frag_atcharge))

            tmp_topo_l = []
            for item in self.Supra_topo_l:
                at1, at2 = item["atoms"]  # 1-based topo labels
                if at1 in global_to_local and at2 in global_to_local:
                    new_item = {
                        "atoms": [global_to_local[at1], global_to_local[at2]],
                        "value": item["value"]
                    }
                    if "mbo" in item:
                        new_item["mbo"] = item["mbo"]
                    tmp_topo_l.append(new_item)

            frag_coord_list.append(frag_coord)
            frag_ele_list.append(frag_ele)
            frag_atcharge_list.append(frag_atcharge)
            frag_c_float_list.append(frag_c_float)
            frag_topo_list.append(tmp_topo_l)
            frag_global_atomidx_list.append(frag_global_atomidx)

        float_charge_sum = float(sum(frag_c_float_list))

        frag_c_int_list = self.adjust_fragment_charge(frag_c_float_list)

        for frag_ele, frag_charge in zip(frag_ele_list, frag_c_int_list):
            frag_nel = sum(element[e] for e in frag_ele) - frag_charge
            frag_nel_list.append(int(frag_nel))

        if sum(frag_nel_list) != self.totnel:
            diff = self.totnel - sum(frag_nel_list)
            raise ValueError(
                f"Fragment electron count mismatch: got {sum(frag_nel_list)}, "
                f"expected {self.totnel}. Charge assignment may be inconsistent by {diff} e-."
            )

        self.Fragments_coordinate_list = frag_coord_list
        self.Fragments_ele_l_list = frag_ele_list
        self.Fragments_atcharge_list = frag_atcharge_list
        self.Fragments_c_list = frag_c_int_list
        self.Fragments_nel_list = frag_nel_list
        self.Fragments_topo_list = frag_topo_list
        self.Fragments_global_atomidx_list = frag_global_atomidx_list

    def adjust_fragment_charge(self, frag_c_float_list):
        charge_l_flo = np.array(frag_c_float_list)
        int_charges = np.floor(charge_l_flo).astype(int)
        residuals = charge_l_flo - int_charges

        n_to_add = int(self.totc - np.sum(int_charges))
        order = np.argsort(-residuals)

        for i in range(n_to_add):
            int_charges[order[i]] += 1

        return list(int_charges)

        # split fragments from Supramolecule: based on strucGraph get the idx of each frag, coord, ele, atcharge all mapped from Supramolecule,for each fragm,  make sure tot_c = tot charge of all fragments, so as nel,


class TS_setting:
    def __init__(self, Rcoordinate, Rele_l, Pcoordinate, Pele_l, c, m):
        self.Rcoordinate = Rcoordinate
        self.Rele_l = Rele_l
        self.Pcoordinate = Pcoordinate
        self.Pele_l = Pele_l
        self.c = c
        self.m = m

        self.Tcoordinate = None
        self.Tele_l = None
        self.features = None


class TSstructure_NEBTS(TS_setting):
    """
    input: reactant/aligned reactants as R; endxyz as P
    => NEBTS (TS guess)
    => g-xtb (E)
    """

    def _get_TSguess(self, comdir):
        os.chdir(comdir)
        write_xyz(self.Rcoordinate, self.Rele_l, "reactant", c=self.c, m=self.m)
        write_xyz(self.Pcoordinate, self.Pele_l, "product", c=self.c, m=self.m)
        rc = CallORCA(comdir, custom_settings={
            "alpb": storage.setting.xtb_setting.alpb,
            "charge": self.c,
            "m": self.m, }).call_nebts()
        if rc == 0:
            self.Tcoordinate, self.Tele_l = read_coord(join(comdir, "orca_NEB-TS_converged.xyz"))
        return rc

    def run(self, rootfolder, debugtrigger):
        random_name = get_random_folder_name()
        comdir = join(rootfolder, random_name)
        create_folder(comdir)
        rc = self._get_TSguess(comdir)
        if rc != 0 and self.barrierless(comdir):
            rc = 0
            self.features = {
                "Esp_gas": 0,
                "Gtrv": 0,
                "Gsolv_corr": 0,
                "G": 0,
                "barrierless": True
            }
            del_dir(not debugtrigger, rootfolder)
            return rc
        elif rc == 0:
            comdir = join(rootfolder, random_name)
            create_folder(comdir)
            self.features = Singlestruc(self.Tcoordinate, self.Tele_l,
                                        self.c, self.m).process_workflow(comdir, opt=False, debugmode=debugtrigger)
            self.features["barrierless"] = False
        else:
            self.features = {
                "Esp_gas": np.inf,
                "G": np.inf,
                "Gtrv": np.inf,
                "Gsolv_corr": np.inf,
                "barrierless": False
            }
        del_dir(not debugtrigger, rootfolder)
        return rc

    @staticmethod
    def barrierless(comdir):
        with open(join(comdir, "orca.out"), "r") as file:
            for line in file.readlines():
                if "No barrier was found." in line:
                    return True
        return False


class TSstructure_ffits(TS_setting):
    """
    input: reactant/aligned reactants as R; endxyz as P
    => FFits (TS guess)
    => ORCA(xtb) (TS geom refinement)
    => g-xtb (E)
    """

    def _get_TSguess(self, comdir):
        # write reactant and prod in comdir
        os.chdir(comdir)
        write_xyz(self.Rcoordinate, self.Rele_l, "reactant", c=self.c, m=self.m)
        write_xyz(self.Pcoordinate, self.Pele_l, "product", c=self.c, m=self.m)
        rc = CallFFiTS(comdir, custom_settings={"charge": self.c, "uhf": self.m - 1})
        if rc == 0:
            self.Tcoordinate, self.Tele_l = read_coord(join(comdir, "optimized.xyz"))
        return rc

    def _refine_TSguess(self, comdir):
        os.chdir(comdir)
        write_xyz(self.Tcoordinate, self.Tele_l, "coord", c=self.c, m=self.m)
        proc = CallORCA(comdir, custom_settings={"charge": self.c, "m": self.m})
        rc = proc.call_optts()
        if rc == 0 and not proc.has_single_imaginary_mode():
            rc = 1
        return rc

    def run(self, rootfolder, debugtrigger=False):
        random_name = get_random_folder_name()
        comdir = join(rootfolder, random_name, "ffits")
        create_folder(comdir)
        rc = self._get_TSguess(comdir)
        if rc == 0:
            comdir = join(rootfolder, random_name, "optts")
            create_folder(comdir)
            rc = self._refine_TSguess(comdir)
            if rc == 0:
                self.Tcoordinate, self.Tele_l = read_coord(join(comdir, "orca.xyz"))
                comdir = join(rootfolder, random_name)
                create_folder(comdir)
                self.features = Singlestruc(self.Tcoordinate, self.Tele_l,
                                            self.c, self.m).process_workflow(comdir, opt=False, debugmode=debugtrigger)
            else:
                self.features = {
                    "Esp_gas": np.inf,
                    "G": np.inf
                }
        else:
            self.features = {
                "Esp_gas": np.inf,
                "G": np.inf
            }
        del_dir(not debugtrigger, rootfolder)
        return rc


class StrucGraph:

    def __init__(self, coordinate, ele_l, topo_l):

        self.coordinate = coordinate
        self.ele_l = ele_l
        self.topo_l = topo_l
        self.graph = self.build_graph()

    def build_graph(self):
        G = nx.Graph()

        # add node
        for i, (ele, coord) in enumerate(zip(self.ele_l, self.coordinate), start=1):  # bc topo_l are index starts 1
            G.add_node(
                i,
                element=ele,
                coord=np.array(coord, dtype=float)
            )

        # add bonds/topology as edges
        # with topology: self.topo_l is func(get_topo_covaradii); otherwise acc. mbo>0.5
        for item in self.topo_l:
            at1, at2 = item["atoms"]

            mbo = item.get("mbo")

            if mbo is not None and mbo < 0.5:
                continue

            attrs = {"distance": item["value"]}

            if mbo is not None:
                attrs["mbo"] = mbo

            G.add_edge(at1, at2, **attrs)
        return G

    def count_fragments(self):
        return nx.number_connected_components(self.graph)

    def get_fragments(self):
        return list(nx.connected_components(self.graph))

    def bond_exist(self, at_l):
        target = set(at_l)
        exists = any(set(item["atoms"]) == target for item in self.topo_l)
        return exists

    def get_atomic_path_of_2atoms(self, at1, at2):
        try:
            # shortest path based on atomic distance
            atomicpath = nx.shortest_path(
                self.graph,
                source=at1,
                target=at2,
                weight="distance"
            )

            # total distance along the path
            path_length = nx.shortest_path_length(
                self.graph,
                source=at1,
                target=at2,
                weight="distance"
            )

            return atomicpath, path_length

        except (nx.NetworkXNoPath, nx.NodeNotFound):
            return None, None


class Allstruc:
    def __init__(self, model):
        self.model = model
        self._next_id = 0
        self._lock = threading.Lock()
        self.data = {}
        self.strucitems_to_sid = {}  # key = (mb_cut,c,m); val = sid
        self.sid_to_strucitems = {}  # key = sid; val = (mb_cut,c,m)

        # record all supra P guess structure
        self.Psupra_to_scanwd = {}  # key = (mb_FULL,c,m); val = scanwd
        self.scanwd_to_Psupra = {}

    def get(self, sid):
        return self.data.get(sid)

    def placeholde_new_sid(self):
        with self._lock:
            sid = self._next_id
            self._next_id += 1
        s = structure_class.Structure(id=sid).struc_info
        return sid, s

    def register_new_sid(self, sid, s, features, movefolder=False):
        self.data[sid] = s
        for k, v in features.items():
            if hasattr(s.general, k):
                setattr(s.general, k, v)
        s.energy.Esp_gas = features["Esp_gas"]
        s.energy.G = features["G"]
        s.energy.Gsolv_corr = features["Gsolv_corr"]
        s.energy.Gtrv = features["Gtrv"]
        mb = features["mb"]
        c = features["charge"]
        m = features["multiplicity"]
        if not "unknow" in mb:
            mbcut = " | ".join(features["mb"].split(" | ")[:storage.setting.searching_setting.general.mbcut_length])
            self.strucitems_to_sid[(mbcut, c, m)] = sid
            self.sid_to_strucitems[sid] = (mbcut, c, m)
        else:
            self.strucitems_to_sid[(mb, c, m)] = sid
            self.sid_to_strucitems[sid] = (mb, c, m)

        if movefolder:
            if not exists(features.get("strucwd")):
                raise ValueError(f"strucwd for {features.get('strucwd')} not exist")
            tmpwd = features["strucwd"]
            tarwd = join(storage.setting.searching_setting.general.strucwd, str(sid))
            if exists(tarwd):
                shutil.rmtree(tarwd)
            shutil.move(tmpwd, tarwd)

        return sid

    @staticmethod
    def update_current_structure_info(features, refid):
        s = storage.allstruc.get(refid)
        for k, v in features.items():
            if hasattr(s.general, k):
                setattr(s.general, k, v)

        try:
            tmpwd = features["strucwd"]
            tarwd = join(storage.setting.searching_setting.general.strucwd, str(refid))
            if exists(tarwd):
                shutil.rmtree(tarwd)
            shutil.move(tmpwd, tarwd)
        except:
            pass

    def register_new_SupraMol(self, mb_full, c, m, scanwd):
        self.Psupra_to_scanwd[(mb_full, c, m)] = scanwd
        self.scanwd_to_Psupra[scanwd] = (mb_full, c, m)


def optimize_single_fragment(coordinate, ele_l, charge, nel, mb, rootfolder, debugmode=False):
    import numpy as np
    from os.path import join

    m_l = [1, 3] if nel % 2 == 0 else [2]
    E_sp = np.inf
    best_features = None
    best_comdir = None
    ideal = False

    for m in m_l:
        random_name = get_random_folder_name()
        comdir = join(rootfolder, random_name)

        tmpfeatures = Singlestruc(
            coordinate, ele_l, charge, m, mb
        ).process_workflow(rootfolder=comdir, ideal=ideal, debugmode=debugmode)
        if tmpfeatures["Esp_gas"] < E_sp:
            if best_comdir is not None:
                del_dir(True, best_comdir)
            E_sp = tmpfeatures["Esp_gas"]
            best_features = tmpfeatures
            best_comdir = comdir
        elif np.isinf(tmpfeatures["Esp_gas"]) and best_features is None:
            best_features = tmpfeatures

        # del_dir(not debugmode, rootfolder)
    return best_features


def process_fragment_recursive(supra_info, rootfolder, depth=0, maxdepth=2, debugmode=False):
    """
        Split first, then optimize each fragment, then recurse until each final
        processed object is a single-fragment structure.

        Returns
        -------
        list[dict]
            A list of processed feature dicts for all terminal single-structure fragments.
        """
    from molbar.barcode import get_molbar_from_coordinates

    coordinate = supra_info["coordinate"]
    ele_l = supra_info["ele_l"]
    charge = supra_info["charge"]
    nel = supra_info["nel"]
    atcharge = supra_info["atcharge"]
    atomic_connection = supra_info["atomic_connection"]
    global_atomidx = supra_info.get("global_atomidx")
    if global_atomidx is None:
        global_atomidx = list(range(1, len(ele_l) + 1))
    mb = supra_info.get("mb")
    if mb is None:
        mb = get_molbar_from_coordinates(coordinate, ele_l)
    # separate
    sup = Supramolecule(
        Supra_coordinate=coordinate,
        Supra_ele_l=ele_l,
        Supra_atcharge=atcharge,
        Supra_topo_l=atomic_connection,
        totnel=nel,
        totc=charge,
        global_atomidx=global_atomidx,
    )
    sup.split()

    if len(sup.Fragments_coordinate_list) > 1 and depth < maxdepth:
        out = []
        for sub_coord, sub_ele, sub_c, sub_nel, sub_topo_l, sub_atcharge, sub_global_atomidx, in zip(
                sup.Fragments_coordinate_list,
                sup.Fragments_ele_l_list,
                sup.Fragments_c_list,
                sup.Fragments_nel_list,
                sup.Fragments_topo_list,
                sup.Fragments_atcharge_list,
                sup.Fragments_global_atomidx_list
        ):
            out.extend(process_fragment_recursive(
                {
                    "coordinate": sub_coord,
                    "ele_l": sub_ele,
                    "charge": sub_c,
                    "nel": sub_nel,
                    "atcharge": sub_atcharge,
                    "atomic_connection": sub_topo_l,
                    "global_atomidx": sub_global_atomidx
                },
                rootfolder,
                depth=depth + 1,
                maxdepth=maxdepth, debugmode=debugmode
            ))
        return out
    # for the single-structure
    opt_info = optimize_single_fragment(coordinate, ele_l, charge, nel, mb, rootfolder, debugmode=debugmode)
    opt_coord = opt_info["coordinate"]
    opt_ele_l = opt_info["ele_l"]
    opt_charge = opt_info["charge"]
    opt_nel = opt_info["nel"]
    opt_atcharge = opt_info["atcharge"]
    opt_topo = opt_info["atomic_connection"]
    opt_info["global_atomidx"] = global_atomidx
    if np.isinf(opt_info["Esp_gas"]):
        # error structure
        return [opt_info]

    post_sup = Supramolecule(
        Supra_coordinate=opt_coord,
        Supra_ele_l=opt_ele_l,
        Supra_atcharge=opt_atcharge,
        Supra_topo_l=opt_topo,
        totnel=opt_nel,
        totc=opt_charge,
        global_atomidx=global_atomidx
    )
    post_sup.split()

    if len(post_sup.Fragments_coordinate_list) > 1 and depth < maxdepth:
        out = []
        for sub_coord, sub_ele, sub_c, sub_nel, sub_topo_l, sub_atcharge, sub_global_atomidx in zip(
                post_sup.Fragments_coordinate_list,
                post_sup.Fragments_ele_l_list,
                post_sup.Fragments_c_list,
                post_sup.Fragments_nel_list,
                post_sup.Fragments_topo_list,
                post_sup.Fragments_atcharge_list,
                post_sup.Fragments_global_atomidx_list
        ):
            out.extend(process_fragment_recursive(
                {
                    "coordinate": sub_coord,
                    "ele_l": sub_ele,
                    "charge": sub_c,
                    "nel": sub_nel,
                    "atcharge": sub_atcharge,
                    "atomic_connection": sub_topo_l,
                    "global_atomidx": sub_global_atomidx
                },
                rootfolder,
                depth=depth + 1,
                maxdepth=maxdepth, debugmode=debugmode
            ))
        return out

    return [opt_info]


def map_fragments_topo_l_back_to_Supramol(products_info):
    map_topo_l = []
    for Pfeatures in products_info:
        local_topo_l = Pfeatures["atomic_connection"]
        global_atomidx = Pfeatures["global_atomidx"]
        for item in local_topo_l:
            local_a, local_b = item["atoms"]

            global_a = global_atomidx[local_a - 1]
            global_b = global_atomidx[local_b - 1]

            new_item = dict(item)
            new_item["atoms"] = sorted([global_a, global_b])

            map_topo_l.append(new_item)
    return map_topo_l


def check_target_bond_change(at_l, map_topo_l, rxn_type):
    bond_set = {
        tuple(sorted(item["atoms"]))
        for item in map_topo_l
    }

    # convert 0-based -> 1-based
    target_bonds = {
        tuple(sorted([a + 1, b + 1]))
        for a, b in at_l
    }

    is_formation = rxn_type in ["interA", "intraA"]

    if is_formation:
        return all(pair in bond_set for pair in target_bonds)
    else:
        return all(pair not in bond_set for pair in target_bonds)
