import os
from os.path import exists, join
import numpy as np
from aRST.script.callsys import CallxTB, CallORCA, CallgXTB, CallMolBar
import aRST.script.initialization as storage
from aRST.script.toolbox import create_folder, read_coord, read_atcharge, write_xyz, read_bo, get_topo
import threading
from aRST.script import structure_class
import networkx as nx
from aRST.script.element import element
import shutil


class Singlestruc:
    # complest electronic properties of a given single structure
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

    def _run_opt(self, comdir, engine="xtb"):
        self._prepare_workdir_and_write_coord(comdir)
        if engine == "orca":
            proc = CallORCA(comdir, custom_settings={"charge": self.c, "m": self.m, "caltyp": "SP OPT"}).call(
                type="opt")
            if proc == 0:
                self.coordinate, self.ele_l = read_coord(join(comdir, "orca.xyz"))
            else:
                print(f'Warning: struc in {comdir} opt at ORCA failed, continue with unopt struc')
        elif engine == "xtb":
            proc = CallxTB(comdir, custom_settings={"charge": self.c, "uhf": self.m - 1}).opt()
            if proc == 0:
                self.coordinate, self.ele_l = read_coord(join(comdir, "xtbopt.xyz"))
            else:
                proc = CallxTB(comdir, custom_settings={"charge": self.c, "uhf": self.m - 1}).sp()
                print(f'Warning: struc in {comdir} opt at xtb failed, try with unopt struc')
                if proc == 0:
                    return self.coordinate, self.ele_l
                else:
                    self.coordinate, self.ele_l
        else:
            raise ValueError(f"Unknown engine for opt: {engine}")
        return self.coordinate, self.ele_l

    def _run_sp(self, comdir, engine="gxtb"):
        self._prepare_workdir_and_write_coord(comdir)
        if engine == "gxtb":
            proc = CallgXTB(comdir, custom_settings={"charge": self.c, "uhf": self.m - 1})
            # write .charg and .uhf file
            with open(".CHRG", "w") as file:
                file.write(f"{self.c}\n")
            with open(".UHF", "w") as file:
                file.write(f"{self.m - 1}\n")
            status = proc.sp()
            if status == 0:
                energy = proc.get_Esp_from_out()
            else:
                energy = np.inf
                print(f'Warning: struc in {comdir} sp at gxtb failed, continue anyhow')

        elif engine == "orca":
            proc = CallORCA(comdir, custom_settings={"charge": self.c, "m": self.m, "caltyp": "SP"})
            staus = proc.call(type="sp")
            if staus == 0:
                energy = proc.get_energy_from_out()
            else:
                energy = np.inf
                print(f'Warning: struc in {comdir} sp at ORCA failed, continue anyhow')
        else:
            raise ValueError(f"Unknown engine for sp: {engine}")
        if energy is None:
            energy = np.inf
            print(f'Warning: struc in {comdir} sp prop fails, structure energy will be set as inf')
        return energy

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

    def process_workflow(self, rootfolder, ideal=True):
        # ideal->opt/sp -> collect res

        from aRST.script.toolbox import get_topo_covaradii
        if len(self.ele_l) > 1 and ideal:
            topo_l_simple = get_topo_covaradii(self.coordinate, self.ele_l)
            G = StrucGraph(self.coordinate, self.ele_l, topo_l_simple).graph
            if nx.is_connected(G):
                self._idealize(topo_l_simple)

        ## opt
        dftopt = storage.setting.read.general.dftopt or storage.setting.orca_setting.dftopt
        if len(self.ele_l) > 1:
            # gfn2 opt
            comdir = join(rootfolder, "xtb")
            self._run_opt(comdir, engine="xtb")
            if dftopt:
                comdir = join(rootfolder, "orca")
                self._run_opt(comdir, engine="orca")
        ## sp
        ## gfn2 sp is necessary
        if dftopt:
            comdir = join(rootfolder, "xtb")
            self._run_sp(comdir, engine="xtb")

        if storage.setting.searching_setting.general.sp_method == "dft":
            comdir = join(rootfolder, "orca")
            if dftopt:
                # read energy from opt proc
                energy = CallORCA(comdir).get_energy_from_out()
            else:
                energy = self._run_sp(comdir, engine="orca")
        elif storage.setting.searching_setting.general.sp_method == "gxtb":
            comdir = join(rootfolder, "gxtb")
            energy = self._run_sp(comdir, engine="gxtb")
        else:
            raise ValueError("sp_method keyword wrong")

        # collect
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
            "Esp": energy,
            "strucwd": rootfolder
        }
        return features

    @staticmethod
    def read_features_from_struc_folder(rootfolder, charge, m):
        # update coord and ele_l: dftopt>gxtb>xtb
        coordinate, ele_l, nel, atcharge, Esp, constrain_l, mb = (None,) * 7

        if exists(join(rootfolder, "orca", "orca.xyz")):
            coordinate, ele_l = read_coord(join(rootfolder, "orca", "orca.xyz"))


        elif exists(join(rootfolder, "xtb", "xtbopt.xyz")):
            coordinate, ele_l = read_coord(join(rootfolder, "xtb", "xtbopt.xyz"))


        elif exists(join(rootfolder, "orca", "coord.xyz")):  # dft sp
            coordinate, ele_l = read_coord(join(rootfolder, "orca", "coord.xyz"))


        elif exists(join(rootfolder, "gxtb", "coord.xyz")):  # gxtb sp
            coordinate, ele_l = read_coord(join(rootfolder, "gxtb", "coord.xyz"))


        else:
            return {}, False

        if storage.setting.searching_setting.general.sp_method == "dft":
            Esp = CallORCA(join(rootfolder, "orca")).get_energy_from_out()
            if m is None:
                m = CallORCA(join(rootfolder, "orca")).get_m()
        else:
            Esp = CallgXTB(join(rootfolder, "gxtb")).get_Esp_from_out()
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
            "atcharge": read_atcharge(join(rootfolder, "xtb", "charges")),
            "mb": mb,
            "atomic_connection": constrain_l,
            "Esp": Esp,
            "strucwd": rootfolder
        }
        if any(v is None for v in features.values()):
            return {}, False

        return features, True


class Supramolecule:
    # analyse the possible Supramolecule
    def __init__(self, Supra_coordinate, Supra_ele_l, Supra_atcharge, Supra_topo_l, totc, totnel):
        self.Supra_coordinate = Supra_coordinate
        self.Supra_ele_l = Supra_ele_l
        self.Supra_atcharge = Supra_atcharge
        self.Supra_topo_l = Supra_topo_l
        self.totc = totc
        self.totnel = totnel
        self.Sgraph = StrucGraph(coordinate=Supra_coordinate, ele_l=Supra_ele_l, topo_l=Supra_topo_l).graph

        self.Fragments_coordinate_list = []
        self.Fragments_ele_l_list = []
        self.Fragments_c_list = []
        self.Fragments_atcharge_list = []
        self.Fragments_nel_list = []
        self.Fragments_atomidx_list = []
        self.Fragments_topo_list = []

    def split(self):
        components = list(nx.connected_components(self.Sgraph))

        frag_coord_list = []
        frag_ele_list = []
        frag_atcharge_list = []
        frag_c_float_list = []
        frag_nel_list = []
        frag_topo_list = []

        # print("sum Supra_atcharge =", float(np.sum(self.Supra_atcharge)))
        # print("totc =", self.totc)
        for comp in components:

            comp_sorted = sorted(comp)
            atom_idxs0 = [i - 1 for i in comp_sorted]

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
        for item in self.topo_l:
            at1, at2 = item["atoms"]
            attrs = {"distance": item["value"]}
            mbo = item.get("mbo")
            if mbo is not None:
                attrs["mbo"] = mbo
            G.add_edge(at1, at2, **attrs)

        return G

    def count_fragments(self):
        return nx.number_connected_components(self.graph)

    def get_fragments(self):
        return list(nx.connected_components(self.graph))


class Allstruc:
    def __init__(self, model):
        self.model = model
        self._next_id = 0
        self._lock = threading.Lock()
        self.data = {}
        self.strucitems_to_sid = {}  # key = (mb_cut,c,m); val = sid
        self.sid_to_strucitems = {}  # key = sid; val = (mb_cut,c,m)

        # record all supra P guess structure
        self.Psupra_to_rxnid = {}  # key = (mb_FULL,c,m); val = rxnid
        self.rxnid_to_Psupra = {}

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
        s.energy.sp = features["Esp"]
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
                raise ValueError(f"strucwd for {sid} not exist")
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

    def register_new_SupraMol(self, mb_full, c, m, refid):
        self.Psupra_to_rxnid[(mb_full, c, m)] = refid
        self.rxnid_to_Psupra[refid] = (mb_full, c, m)

def optimize_fragment(coordinate, ele_l, charge, nel, mb, rootfolder):
    from aRST.script.toolbox import get_random_folder_name
    import numpy as np
    from os.path import join

    m_l = [1, 3] if nel % 2 == 0 else [2]
    E_sp = np.inf
    best_features = None
    ideal = False

    for m in m_l:
        random_name = get_random_folder_name()
        comdir = join(rootfolder, random_name)

        tmpfeatures = Singlestruc(
            coordinate, ele_l, charge, m, mb
        ).process_workflow(rootfolder=comdir, ideal=ideal)

        if tmpfeatures["Esp"] < E_sp:
            E_sp = tmpfeatures["Esp"]
            best_features = tmpfeatures

    return best_features
def process_fragment_recursive(supra_info, rootfolder, depth=0,maxdepth=2):
    """
        Split first, then optimize each fragment, then recurse until each final
        processed object is a single-fragment structure.

        Returns
        -------
        list[dict]
            A list of processed feature dicts for all terminal single-structure fragments.
        """
    from aRST.script.toolbox import get_random_folder_name
    from molbar.barcode import get_molbar_from_coordinates

    coordinate = supra_info["coordinate"]
    ele_l = supra_info["ele_l"]
    charge = supra_info["charge"]
    nel = supra_info["nel"]
    atcharge = supra_info["atcharge"]
    atomic_connection = supra_info["atomic_connection"]
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
        totc=charge
    )
    sup.split()


    if depth >= maxdepth:
        return [optimize_fragment(coordinate, ele_l, charge, nel, mb, rootfolder)]

    if len(sup.Fragments_coordinate_list) == 1:
        return [optimize_fragment(coordinate, ele_l, charge, nel, mb, rootfolder)]

    out = []
    for sub_coord, sub_ele, sub_c, sub_nel,sub_topo_l,sub_atcharge in zip(
            sup.Fragments_coordinate_list,
            sup.Fragments_ele_l_list,
            sup.Fragments_c_list,
            sup.Fragments_nel_list,
            sup.Fragments_topo_list,
            sup.Fragments_atcharge_list,
    ):
        frag_info = {
            "coordinate":sub_coord,
            "ele_l":sub_ele,
            "charge":sub_c,
            "nel":sub_nel,
            "atcharge":sub_atcharge,
            "atomic_connection":sub_topo_l
        }

        out.extend(
            process_fragment_recursive(
                frag_info,
                rootfolder,
            )
        )

    return out