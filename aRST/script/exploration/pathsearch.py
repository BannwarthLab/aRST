from aRST.script.toolbox import create_folder, write_xyz, read_coord, read_bo, read_atcharge, get_topo,add_noisy_to_coord
from aRST.script.callsys import CallxTB
import aRST.script.initialization as storage
from os.path import join
from molbar.barcode import get_molbar_from_coordinates
from aRST.script.restart.checkpoint import checkpoint_scan_search,checkpoint_mbff_search
import numpy as np
import shutil


class SinglePathSearch:
    def __init__(self, coordinate_init, ele_l_init, u, charge, comdir, nat_1=None, reactants=None):
        self.coordinate_init = coordinate_init
        self.ele_l_init = ele_l_init
        self.comdir = comdir
        self.u = u
        self.charge = charge
        self.nat_1 = nat_1  # placeholder for inter mol, the at number of the first reactant
        self.reactants = reactants

    def _short_md_for_product_struc_guess(self, coordinate, ele_l, constrained_bonds, charge, uhf):
        mdwd = join(self.comdir, "md")
        create_folder(mdwd)
        write_xyz(coordinate, ele_l, "coord")
        mdtemp = storage.setting.xtb_setting.scan_etemp
        mdsetting = storage.setting.searching_setting.md
        xtbmdparse = CallxTB(mdwd, custom_settings={"charge": charge, "uhf": uhf,
                                                    "mdetemp": mdtemp,
                                                    "mdtime": mdsetting.time,
                                                    "mdstep": mdsetting.step,
                                                    "mdshake": mdsetting.shake,
                                                    "mddump": mdsetting.dump,
                                                    "etemp": storage.setting.xtb_setting.scan_etemp
                                                    # sp after md
                                                    })
        status = xtbmdparse.short_md(constrained_bonds=constrained_bonds)
        if status != 0:
            raise ValueError(f"md has error, check {mdwd}")
        # get energy min struc in the md simulation
        # coord can't be single at in this case
        coordinate, ele_l, emd = xtbmdparse.get_Emin_geom_from_mdtrj(nat=len(ele_l))
        write_xyz(coordinate, ele_l, "emincoord")
        # doing a sp in the same folder for sp information
        status = xtbmdparse.sp()
        if status != 0:
            raise ValueError(f"sp after md has error, check {mdwd}")
        bo = read_bo(join(mdwd, "wbo"))
        constrain_l = get_topo(coordinate, ele_l, bo)
        atcharge = read_atcharge(join(mdwd, "charges"))
        return bo, constrain_l, atcharge, coordinate, ele_l

    def _add_constrain_for_loose_weakest_bond(self, constrain_setting, generalsetting):

        # searching weakest bond
        from aRST.script.spr import SDP
        from aRST.script.geom import Coord

        # searching weakest bond
        allbonds = SDP(self.reactants,
                       exclude_elements=generalsetting.exclude_elements,
                       exclude_pairs=generalsetting.exclude_pairs
                       ).allbonds
        allbonds = sorted(allbonds, key=lambda x: x["mbo"])
        if not allbonds:
            return constrain_setting

        weakest_bond_entry = allbonds[0]
        if weakest_bond_entry["mbo"] > 0.9:
            return constrain_setting

        # decide it is inra or inter reactants
        if len(self.reactants) == 1:
            at1, at2 = weakest_bond_entry["atl"]
            at1, at2 = int(at1), int(at2)
        elif len(self.reactants) == 2 and self.nat_1 is not None:
            at1, at2 = weakest_bond_entry["atl"]
            at1, at2 = int(at1) + self.nat_1, int(at2) + self.nat_1
        else:
            raise ValueError("func: _add_constrain_for_loose_weakest_bond, len(reactants) doesnt match 1 or 2")
        dd = np.round(Coord(self.coordinate_init, self.ele_l_init).get_ats_distance(int(at1), int(at2)), 3)
        constrain_setting["constrain"]["loose_bond"] = at1 + 1, at2 + 1, dd
        return constrain_setting

    @checkpoint_scan_search(dirname="checkpoint/explore")
    def run_one_scan_search(self, constrain_setting, rxn_type, reactantinfo, nel):
        searchingsetting = storage.setting.searching_setting
        generalsetting = storage.setting.searching_setting.general
        if searchingsetting.inter.loose_bond or searchingsetting.intra.loose_bond:
            constrain_setting = self._add_constrain_for_loose_weakest_bond(constrain_setting, generalsetting)

        create_folder(self.comdir)
        # write coord
        write_xyz(self.coordinate_init, self.ele_l_init, "coord")

        # call xtb scan
        status = CallxTB(self.comdir,
                         custom_settings={"charge": self.charge, "uhf": self.u,
                                          "etemp": storage.setting.xtb_setting.scan_etemp}).scan(constrain_setting)
        # receive results
        results = None
        if status == 0:
            endxyz, ele_l_end = read_coord(join(self.comdir, "xtbopt.xyz"))
            bo = read_bo(join(self.comdir, "wbo"))
            constrain_l = get_topo(endxyz, ele_l_end, bo)
            atcharge = read_atcharge(join(self.comdir, "charges"))
            if generalsetting.sep_mode != 0:
                # do a unconstarined opt
                optwd = join(self.comdir, "opt")
                create_folder(optwd)
                write_xyz(endxyz, ele_l_end, "coord")
                status2 = CallxTB(optwd, custom_settings={"charge": self.charge, "uhf": self.u}).opt()
                if status2 == 0:
                    endxyz, ele_l_end = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l = get_topo(endxyz, ele_l_end, bo)
                    atcharge = read_atcharge(join(optwd, "charges"))

            # do relax opt optional 
            if generalsetting.sep_mode != 0:
                # do a unconstarined opt
                optwd = join(self.comdir, "opt")
                create_folder(optwd)
                shutil.copyfile(join(self.comdir, 'xtbopt.xyz'), join(optwd, 'coord.xyz'))
                status2 = CallxTB(optwd, custom_settings={"charge": self.charge, "uhf": self.u}).opt()
                if status2 == 0:
                    endxyz, ele_l_end = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l = get_topo(endxyz, ele_l_end, bo)
                    atcharge = read_atcharge(join(optwd, "charges"))

            # md optional
            if generalsetting.md:
                push_l = constrain_setting["scan"].get("distance")
                bo, constrain_l, atcharge, endxyz, ele_l_end = self._short_md_for_product_struc_guess(endxyz,
                                                                                                      ele_l_end, push_l,
                                                                                                      self.charge,
                                                                                                      self.u)

            endxyz_mb = get_molbar_from_coordinates(endxyz, ele_l_end)
            results = [
                reactantinfo,
                {"atcharge": atcharge,
                 "charge": self.charge,
                 "uhf": self.u,
                 "nel": nel,
                 "coordinate": endxyz,
                 "ele_l": ele_l_end,
                 "bo": bo,
                 "atomic_connection": constrain_l,
                 "mb": endxyz_mb},  # endxyz info
                {"type": rxn_type,
                 "method": "scan",
                 "push_pair":constrain_setting["scan"]["distance"]
                 }  # reaction info
            ]

        return results,status

    @checkpoint_mbff_search(dirname="checkpoint/explore")
    def run_xtbopt_for_mbff_search(self, atom_connection, push_l, rxn_type, reactantinfo, nel):

        searchingsetting = storage.setting.searching_setting
        generalsetting = storage.setting.searching_setting.general

        create_folder(self.comdir)
        coordinate_noisy = add_noisy_to_coord(self.coordinate_init, scale=1e-4)
        write_xyz(coordinate_noisy, self.ele_l_init, "coord_mb")
        # constrained opt in gfn2
        newbond = [{"atoms":(item[0],item[1]),"value":item[-1]} for item in push_l]
        status = CallxTB(self.comdir,
                         custom_settings={"charge": self.charge, "uhf": self.u,
                                          "k": storage.setting.xtb_setting.k,
                                          "coordfn": "coord_mb",
                                          "etemp": storage.setting.xtb_setting.scan_etemp
                                          }).constrained_opt(
            newbond=newbond)
        results = None
        if status == 0:
            endxyz, ele_l_end = read_coord(join(self.comdir, "xtbopt.xyz"))
            bo = read_bo(join(self.comdir, "wbo"))
            constrain_l = get_topo(endxyz, ele_l_end, bo)
            atcharge = read_atcharge(join(self.comdir, "charges"))

            if generalsetting.sep_mode != 0:
                # do a unconstarined opt
                optwd = join(self.comdir, "opt")
                create_folder(optwd)
                shutil.copyfile(join(self.comdir, 'xtbopt.xyz'), join(optwd, 'coord.xyz'))
                status2 = CallxTB(optwd, custom_settings={"charge": self.charge, "uhf": self.u}).opt()
                if status2 == 0:
                    endxyz, ele_l_end = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l = get_topo(endxyz, ele_l_end, bo)
                    atcharge = read_atcharge(join(optwd, "charges"))
            if generalsetting.md:
                bo, constrain_l, atcharge, endxyz, ele_l_end = self._short_md_for_product_struc_guess(endxyz,
                                                                                                      ele_l_end, push_l,
                                                                                                      self.charge,
                                                                                                      self.u)

            endxyz_mb = get_molbar_from_coordinates(endxyz, ele_l_end)
            results = [
                reactantinfo,
                {"atcharge": atcharge,
                 "charge": self.charge,
                 "uhf": self.u,
                 "nel":nel,
                 "coordinate": endxyz,
                 "ele_l": ele_l_end,
                 "bo": bo,
                 "atomic_connection": constrain_l,
                 "mb": endxyz_mb},  # endxyz info
                {"type": rxn_type,
                 "method": "mbff",
                 "push_pair": push_l
                 }  # reaction info
            ]
        return results,status
