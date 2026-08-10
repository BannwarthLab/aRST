import sys
import traceback
from os.path import join, exists
import os
import shutil
import subprocess
import re
from aRST.script.toolbox import mydict, read_coord
import numpy as np


def remove_file(old_path, new_path):
    filelist = os.listdir(old_path)

    for file in filelist:
        src = join(old_path, file)
        dst = join(new_path, file)
        shutil.move(src, dst)


class CallORCA:
    @classmethod
    def default_settings(cls, orcasetting,mpsetting):
        """Dynamically generates default settings based on the `setting` class."""
        return {
            "functional": orcasetting.functional,
            "basis": orcasetting.basis,
            "charge": 0,
            "m": 1,
            "fn": "coord",
            "caltyp": "SP",
            "caltyp2": "RKS",  # This can be updated dynamically if needed
            "cpcm": orcasetting.cpcm,  # only use cpcm when opt
            "ncpu": mpsetting.orca_nprocs,
            "maxcore": orcasetting.maxcore,
        }

    def __init__(self, wd, command_d=None, custom_settings=None):

        self.wd = wd
        if command_d is not None:
            self.cmd = mydict(command_d)
        else:
            import aRST.script.initialization as storage
            command_d = self.default_settings(storage.setting.orca_setting,storage.setting.mp_setting)
            if custom_settings is not None:
                command_d.update(custom_settings)

            command_d["caltyp2"] = "RKS" if command_d["m"] == 1 else "UKS"
            self.cmd = mydict(command_d)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        self.ncpu = self.cmd.get("ncpu", 1)

    def orcasubprocess(self, cmd_string, writeenergy=True, removegbw=True):
        out_path = f'{self.wd}/orca.out'
        err_path = f'{self.wd}/orca_err.out'

        with open(out_path, 'w') as stdout_file, open(err_path, 'w') as stderr_file:
            p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        if rc == 0:
            if os.path.exists(err_path) and os.path.getsize(err_path) == 0:
                os.remove(err_path)
            if writeenergy:
                with open(out_path, 'r') as f:
                    lines = f.readlines()
                    for line in reversed(lines):
                        if 'FINAL SINGLE' in line:
                            with open(f'{self.wd}/ENERGY', 'w') as f2:
                                f2.write(str(re.findall(r'\-?\d+\.?\d*', line)[0]))
                            break
        extensions_to_remove = ['.densities', '.opt']
        keep_list = ['.xyz', '.inp', '.in', '.out', '.input', '.molden', '.gbw']
        if removegbw:
            keep_list.remove('.gbw')
        for filename in os.listdir(f'{self.wd}'):
            if filename == "ENERGY":
                continue
            if not any(filename.endswith(ext) for ext in keep_list):
                os.remove(filename)

        return rc

    def generate_cmd_string(self, type="sp"):
        command_string = "! "
        command_string += f"{self.cmd.get_value('functional', '')} "
        if self.cmd.get_value('basis', ''):
            command_string += f"{self.cmd.get_value('basis', '')} "
        command_string += f"{self.cmd.get_value('caltyp', '')} "
        command_string += f"{self.cmd.get_value('caltyp2', '')} "
        if type == "sp":
            command_string += f"\n"
        elif self.cmd.get_value('cpcm', ''):
            command_string += f"CPCM({self.cmd.get_value('cpcm', '')}) \n"

        if isinstance(self.ncpu, int) and self.ncpu > 1:
            if self.cmd.get_value('ncpu', 1):
                command_string += f"%PAL NPROCS {int(self.ncpu)} END \n"
                command_string += f"%maxcore {self.cmd.get_value('maxcore', 1000)} \n"

        c = self.cmd.get_value('charge')
        m = self.cmd.get_value('m')
        fn = self.cmd.get_value('fn')
        if all(v is not None for v in [fn, c, m]):
            command_string += f"* xyzfile {c} {m} {fn}.xyz "
        else:
            raise Exception(f"{fn}.xyz in {self.wd} has wrong c,m,fn parameter with {[fn, c, m]}")
        return command_string

    def call(self, removegbw=True, type="sp"):
        command_string = self.generate_cmd_string(type=type)
        with open(f'{self.wd}/orca.inp', 'w') as f:
            f.write(command_string)

        cmd = f"{os.environ['ORCA_BIN']} orca.inp"
        rc = self.orcasubprocess(cmd, removegbw=removegbw)
        return rc

    def get_geomres(self):
        if exists(join(self.wd, 'orca.xyz')):
            coordinate, ele_l = read_coord(join(self.wd, 'orca.xyz'))
        else:
            coordinate, ele_l = read_coord(join(self.wd, 'coord.xyz'))
        return coordinate, ele_l

    def get_energy_from_out(self, out="orca.out"):
        Esp = np.inf
        if exists(join(self.wd, "ENERGY")):
            with open(join(self.wd, "ENERGY"), 'r') as file:
                Esp = float(file.readlines()[0])
        elif exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                lines = file.readlines()
            for ln in reversed(lines):
                if "FINAL SINGLE POINT ENERGY" in ln:
                    Esp = float(ln.split()[-1])
                    break
        return Esp

    def get_Gtrv_from_out(self, out="orca.out"):
        G_trv = np.inf
        if exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                lines = file.readlines()
            for ln in reversed(lines):
                if "G-E(el)" in ln:
                    G_trv = float(ln.split()[-4])
                    return G_trv
        return G_trv

    def get_ZPE_from_out(self, out="orca.out"):
        ZPE= np.inf
        if exists(join(self.wd, out)):
            with open(join(self.wd, out), "r") as file:
                lines = file.readlines()
            for ln in reversed(lines):
                if "Zero point energy" in ln:
                    parts = ln.split()
                    ZPE = float(parts[parts.index("Eh") - 1])
                    return ZPE
        return ZPE

    def get_m(self, file="orca.inp"):
        m = None
        if exists(join(self.wd, file)):
            with open(join(self.wd, file), 'r') as file:
                lines = file.readlines()
            for ln in lines:
                if "xyz" in ln:
                    parts = ln.split()
                    m = int(parts[3])
        return m

    def call_nebts(self):
        # write input
        with open(f'{self.wd}/orca.inp', 'w') as f:
            f.write(f"""
!XTB NEB-TS ALPB({self.cmd.get_value('alpb')})
%NEB
 NEB_END_XYZFILE "product.xyz"
END
%PAL NPROCS {self.cmd.get_value('ncpu')} END
%maxcore {self.cmd.get_value('maxcore')}
*XYZFILE {self.cmd.get_value('charge')} {self.cmd.get_value('m')} reactant.xyz           
            """)

        cmd = f"{os.environ['ORCA_BIN']} orca.inp"

        rc = self.orcasubprocess(cmd, removegbw=False)
        # check converge coord.xyz
        if not os.path.exists(join(self.wd, "orca_NEB-TS_converged.xyz")):
            rc = 1
        return rc

    def call_trv(self):
        with open(f'{self.wd}/orca.inp', 'w') as f:
            f.write(f"""
!XTB FREQ ALPB({self.cmd.get_value('alpb')})
%PAL NPROCS {self.cmd.get_value('ncpu')} END
%maxcore {self.cmd.get_value('maxcore')}
*XYZFILE {self.cmd.get_value('charge')} {self.cmd.get_value('m')} {self.cmd.get_value('fn')}.xyz
""")
        cmd = f"{os.environ['ORCA_BIN']} orca.inp"

        rc = self.orcasubprocess(cmd, removegbw=False)

        return rc
    def call_optts(self):
        with open(f'{self.wd}/orca.inp', 'w') as f:
            f.write(f"""
!XTB optts FREQ ALPB({self.cmd.get_value('alpb')})
%PAL NPROCS {self.cmd.get_value('ncpu')} END
%maxcore {self.cmd.get_value('maxcore')}
*XYZFILE {self.cmd.get_value('charge')} {self.cmd.get_value('m')} coord.xyz
""")
        cmd = f"{os.environ['ORCA_BIN']} orca.inp"

        rc = self.orcasubprocess(cmd, removegbw=False)

        return rc
    @staticmethod
    def has_single_imaginary_mode(out="orca.out"):
        with open(out, "r") as f:
            text = f.read()

        # isolate vibrational frequency block
        match = re.search(
            r"VIBRATIONAL FREQUENCIES(.*?)(?:NORMAL MODES|THERMOCHEMISTRY|$)",
            text,
            re.S,
        )

        if match is None:
            return False

        vib_block = match.group(1)

        n_imag = vib_block.count("***imaginary mode***")

        return n_imag == 1




class CallCrest:
    @classmethod
    def default_settings(cls, xtbsetting):
        """Dynamically generates default settings based on the `setting` class."""
        return {"method": xtbsetting.method,
                "coordfn": "coord",
                "optlevel": xtbsetting.optlevel,
                "charge": 0,
                "uhf": 0,
                "alpb": xtbsetting.alpb,
                "etemp": xtbsetting.struc_etemp,
                }

    def __init__(self, wd, command_d=None, custom_settings=None):
        self.wd = wd
        if command_d is not None:
            self.cmd = mydict(command_d)
        else:
            import aRST.script.initialization as storage
            command_d = self.default_settings(storage.setting.xtb_setting)
            if custom_settings is not None:
                command_d.update(custom_settings)

            self.cmd = mydict(command_d)

    def crestsubprocess(self, cmd_string):
        out_path = f'{self.wd}/crest.out'
        err_path = f'{self.wd}/crest_err.out'
        with open(out_path, 'w') as stdout_file, open(err_path, 'w') as stderr_file:
            p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        if rc == 0 and os.path.exists(err_path) and os.path.getsize(err_path) == 0:
            os.remove(err_path)
        return rc

    def generate_cmd_string(self):
        cmd_string = os.environ['CREST_BIN'] + " "
        coordfn = self.cmd.get_value("coordfn")
        if coordfn is None:
            raise Exception(f"coord in {self.wd} can't find right coordfn name")
        cmd_string += f"{coordfn}.xyz "

        cmd_string += f"--{self.cmd.get_value('method')} " \
            if self.cmd.get_value('method') else ""
        cmd_string += f"-c {self.cmd.get_value('charge')} " \
            if self.cmd.get_value('charge') else ""
        cmd_string += f"--uhf {self.cmd.get_value('uhf')} " \
            if self.cmd.get_value('uhf') else ""
        cmd_string += f"--etemp {self.cmd.get_value('etemp')} " \
            if self.cmd.get_value('etemp') else ""
        cmd_string += f"--alpb {self.cmd.get_value('alpb')} " \
            if self.cmd.get_value('alpb') else ""
        cmd_string += f"{self.cmd.get_value('extracmd')} " \
            if self.cmd.get_value('extracmd') else ""
        return cmd_string

    def depotonation(self):
        cmd_string = self.generate_cmd_string()
        cmd_string += " --deprotonate "
        rc = self.crestsubprocess(cmd_string)
        return rc

    def potonation(self):
        cmd_string = self.generate_cmd_string()
        cmd_string += " --protonate "
        rc = self.crestsubprocess(cmd_string)
        return rc

    def read_update_coord(self, job):
        import numpy as np
        outwd = f'{self.wd}/crest.out'
        with open(outwd, "r") as file:
            lines = file.readlines()
        trigger = False
        startline, end = 0, 0
        for line in reversed(lines):
            if "structure    ΔE(kcal/mol)   Etot(Eh)" in line:
                startline = lines.index(line)
                trigger = True
            if trigger and "Wall Time Summary" in line:
                end = lines.index(line) - 3
                trigger = False
                break
        dE_l = []
        for line in lines[startline + 1:end]:
            nstruc, dE, Etot = line.split()
            dE_l.append(float(dE))
        if min(dE_l) > 0:
            coordinate, ele_l = None, None
        else:
            nstruc = dE_l.index(min(dE_l))
            ensemblewd = f'{self.wd}/deprotonated.xyz' if job == "depotonation" else f'{self.wd}/protonated.xyz'
            coordinate_l = []
            ele_ls = []
            with open(ensemblewd, 'r') as f:
                while True:
                    # Read frame header
                    line = f.readline()
                    if not line:  # End of file
                        break
                    num_atoms = int(line.strip())
                    # Read comment line (contains energy)
                    comment_line = f.readline().strip()

                    elems = np.empty(num_atoms, dtype=str)
                    coords = np.empty((num_atoms, 3), dtype=np.float64)
                    for i in range(num_atoms):
                        parts = f.readline().split()
                        elems[i] = parts[0]
                        coords[i, 0] = float(parts[1])
                        coords[i, 1] = float(parts[2])
                        coords[i, 2] = float(parts[3])

                    ele_ls.append(elems)
                    coordinate_l.append(coords)

            coordinate = coordinate_l[nstruc]
            ele_l = ele_ls[nstruc]
        return coordinate, ele_l


class CallxTB:
    @classmethod
    def default_settings(cls, xtbsetting,mpsetting):
        """Dynamically generates default settings based on the `setting` class."""
        return {"method": xtbsetting.method,
                "coordfn": "coord",
                "optlevel": xtbsetting.optlevel,
                "charge": 0,
                "uhf": 0,
                "alpb": xtbsetting.alpb,
                "etemp": xtbsetting.struc_etemp,
                "nthreads":mpsetting.xtb_threads
                }

    def __init__(self, wd, command_d=None, custom_settings=None, mlversion=False):
        self.wd = wd
        if command_d is not None:
            command_d.update(custom_settings)
            self.cmd = mydict(command_d)

        else:
            import aRST.script.initialization as storage
            command_d = self.default_settings(storage.setting.xtb_setting,storage.setting.mp_setting)
            if custom_settings is not None:
                command_d.update(custom_settings)

            self.cmd = mydict(command_d)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] =str(self.cmd["nthreads"])
        base = f"{self.cmd.get_value('coordfn')}.xyz --{self.cmd.get_value('method')}  -c {self.cmd.get_value('charge')} --uhf {self.cmd.get_value('uhf')} --etemp {self.cmd.get_value('etemp')}"
        if not mlversion:
            self.cmd_string_base = f"{os.environ['XTB_BIN']} {base}  -P {env['OMP_NUM_THREADS']} "
        else:
            self.cmd_string_base = f"{os.environ['XTBML_BIN']} {base} -P {env['OMP_NUM_THREADS']} "

    def xtbsubprocess(self, cmd_string, writeenergy=False):
        out_path = f'{self.wd}/xtb.out'
        err_path = f'{self.wd}/xtb_err.out'
        with open(out_path, 'w') as stdout_file, open(err_path, 'w') as stderr_file:
            p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        if rc == 0:
            if os.path.exists(err_path) and os.path.getsize(err_path) == 0:
                os.remove(err_path)
            if writeenergy:
                with open(err_path, 'r') as f:
                    lines = f.readlines()
                    for line in reversed(lines):
                        if 'TOTAL ENERGY' in line:
                            with open(f'{self.wd}/ENERGY', 'w') as f2:
                                f2.write(str(re.findall(r'\-?\d+\.?\d*', line)[0]))
                            break
        return rc


    def sp_gas(self):
        cmd_string = f"{self.cmd_string_base}"
        rc = self.xtbsubprocess(cmd_string)
        return rc

    def sp_alpb(self):
        cmd_string = f"{self.cmd_string_base}"

        cmd_string +=f"--alpb {self.cmd.get_value('alpb')} "
        rc = self.xtbsubprocess(cmd_string)
        return rc

    def opt(self):
        cmd_string = f"{self.cmd_string_base}--opt "
        cmd_string += f"{self.cmd.get_value('optlevel')} " if self.cmd.get_value('optlevel') else ""
        cmd_string += f"--alpb {self.cmd.get_value('alpb')} "
        rc = self.xtbsubprocess(cmd_string)

        # fallback: retry with high etemp
        if rc != 0:
            oldetmp = self.cmd.get_value('etemp')
            cmd_string_highetemp = f"{self.cmd_string_base}--opt ".replace(f"--etemp {oldetmp}", "--etemp 4000")
            rc = self.xtbsubprocess(cmd_string_highetemp)
            if rc == 0:
                os.rename(f"{self.wd}/xtbopt.xyz",f"{self.wd}/coord.xyz")
                rc = self.xtbsubprocess(cmd_string)
        return rc

    def afo(self):
        cmd_string = f"{self.cmd_string_base}--ml_feature --alpb {self.cmd.get_value('alpb')} "
        p = subprocess.Popen(cmd_string, stdout=open(f'{self.wd}/xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc

    def molden(self):
        cmd_string = f"{self.cmd_string_base}--molden --alpb {self.cmd.get_value('alpb')} "
        p = subprocess.Popen(cmd_string, stdout=open(f'{self.wd}/xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc

    def solv_correct(self):

        ## gas
        cmd_string_gas = f"{self.cmd_string_base}"
        cmd_string_gas +=f"--grad "
        p = subprocess.Popen(cmd_string_gas, stdout=open(f'{self.wd}/gas.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc_gas = p.returncode

        ## solv
        cmd_string_solv = f"{self.cmd_string_base}"
        cmd_string_solv += f"--alpb {self.cmd.get_value('alpb')} bar1M "
        cmd_string_solv += f"--grad "
        p = subprocess.Popen(cmd_string_solv, stdout=open(f'{self.wd}/solv.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc_solv = p.returncode
        return rc_gas,rc_solv


    def scan(self, setting_d):
        # write inp
        with open(f'{self.wd}/scan.inp', 'w') as f:
            f.write(f'$constrain\n'
                    f'force constant={setting_d["constrain"]["force_constant"]}\n')
            if "loose_bond" in setting_d["constrain"].keys():
                at1, at2, dd = setting_d["constrain"]["loose_bond"]
                f.write('distance: %d,%d,%.2f\n' % (int(at1), int(at2), dd * 1.2))

            f.write('$scan\n')
            if setting_d["scan"]["mode"] == "concerted":
                f.write('mode=concerted\n')

            atl = setting_d["scan"]["distance"]
            for item in atl:
                at1, at2, d0, d_end = item
                f.write('distance: %d,%d,%.2f;%.2f,%.2f,%d\n'
                        % (int(at1), int(at2), d0, d0, d_end, setting_d["scan"]["step"]))
            f.write("$end")

        cmd_string = f"{self.cmd_string_base}--opt "
        cmd_string += f"{self.cmd.get_value('optlevel')} " if self.cmd.get_value('optlevel') else ""
        cmd_string += f"--alpb {self.cmd.get_value('alpb')} "
        cmd_string += f"--input scan.inp "

        p = subprocess.Popen(cmd_string, stdout=open(f'{self.wd}/xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        if rc != 0 and (os.path.exists(f'{self.wd}/xtbscan.log') and os.path.exists(f'{self.wd}/xtbopt.xyz')):
            rc = 0
        return rc

    def constrained_opt(self, newbond):
        with open(f'{self.wd}/scan.inp', 'w') as f:
            f.write(f'$constrain\n'
                    f'force constant={self.cmd.get_value("k", "1")}\n')
            for item in newbond:
                at1, at2 = item["atoms"]
                value = item["value"]
                f.write('distance: %d,%d,%.2f\n'
                        % (at1, at2, value))
            f.write("$end")

        cmd_string = f"{self.cmd_string_base}--opt "
        cmd_string += f"{self.cmd.get_value('optlevel')} " if self.cmd.get_value('optlevel') else ""
        cmd_string += f"--alpb {self.cmd.get_value('alpb')} "
        cmd_string += f"--input scan.inp "
        p = subprocess.Popen(cmd_string, stdout=open(f'{self.wd}/xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        if rc != 0 and os.path.exists(f'{self.wd}/xtbopt.xyz'):
            rc = 0
        return rc


    def short_md(self, constrained_bonds=False):
        # constrained_bonds = push_l in the scanning
        with open(f'{self.wd}/md.inp', 'w') as f:
            if constrained_bonds:
                f.write("$constrain\n")
                for item in constrained_bonds:
                    at1 = int(item[0])
                    at2 = int(item[1])
                    d = item[3]
                    f.write(f"distance: {at1}, {at2}, {d:.2f}\n")

            f.write("$md\n"
                    f"time={self.cmd.get_value('mdtime', 10)}\n"  # simulation 10 ps 
                    f"step={self.cmd.get_value('mdstep', 1)}\n"  # every 1 fs for propagation
                    f"temp={self.cmd.get_value('mdetemp', 1000)}\n"  # use set relaxed etemp
                    f"shake={self.cmd.get_value('mdshake', 1)}\n"  # not constring bonds
                    f"dump={self.cmd.get_value('mddump', 50)}"  # save geom every 50 fs
                    )
            # xtb coord.xyz --input md.inp --md

        cmd_string = f"{self.cmd_string_base}--input md.inp --md "
        cmd_string += f"--alpb {self.cmd.get_value('alpb')} "
        cmd_string += f"--input scan.inp "
        p = subprocess.Popen(cmd_string, stdout=open(f'{self.wd}/xtbmd.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc

    def get_Emin_geom_from_mdtrj(self, nat):
        import numpy as np

        # cut geomblock in md.trj
        coordinate_l = []
        ele_ls = []
        energies = []
        trjwd = f'{self.wd}/xtb.trj'

        def parse_xtb_trj(trjwd, nat):
            with open(trjwd, 'r') as f:
                while True:
                    # Read frame header
                    line = f.readline()
                    if not line:  # End of file
                        break
                    num_atoms = int(line.strip())
                    assert num_atoms == nat, f"Unexpected number of atoms: {num_atoms} (expected {nat})"

                    # Read comment line (contains energy)
                    comment_line = f.readline().strip()
                    energy = float(comment_line.split()[1])  # Assumes format: "energy: <value> ..."
                    energies.append(energy)

                    elems = np.empty(num_atoms, dtype=str)
                    coords = np.empty((num_atoms, 3), dtype=np.float64)
                    for i in range(num_atoms):
                        parts = f.readline().split()
                        elems[i] = parts[0]
                        coords[i, 0] = float(parts[1])
                        coords[i, 1] = float(parts[2])
                        coords[i, 2] = float(parts[3])

                    ele_ls.append(elems)
                    coordinate_l.append(coords)
            return coordinate_l, ele_ls, energies

        def filter_outliers(coordinate_l, ele_ls, energies, z_threshold=2.5):
            from scipy.stats import zscore
            z_scores = np.abs(zscore(energies))
            non_outliers = np.where(z_scores < z_threshold)[0]

            clean_coordinate_l = [coordinate_l[i] for i in non_outliers]
            clean_ele_l = [ele_ls[i] for i in non_outliers]
            clean_energies = [energies[i] for i in non_outliers]

            return clean_coordinate_l, clean_ele_l, clean_energies

        coordinate_l, ele_ls, energies = parse_xtb_trj(trjwd, nat)
        clean_coord_l, clean_ele_l, clean_energies = filter_outliers(coordinate_l, ele_ls, energies)
        lowest_idx = np.argmin(clean_energies)
        coordinate_emin = clean_coord_l[lowest_idx]
        ele_l_emin = clean_ele_l[lowest_idx]
        emin = clean_energies[lowest_idx]
        return coordinate_emin, ele_l_emin, emin


    def get_charge_from_out(self, out="xtb.out"):
        c = None
        if exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                lines = file.readlines()

            for ln in lines:
                if "net charge" in ln:
                    c = int(ln.split()[-2])
                    break
        return c

    def get_Esp_from_out(self, out="xtb.out"):
        # unit in eh
        Esp = np.inf

        if exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                lines = file.readlines()

            for ln in reversed(lines):
                if 'TOTAL ENERGY' in ln:
                    Esp = float(ln.split()[-3])
                    break
        return Esp


    def get_Gsolv_corr_from_out(self,solvout="solv.out",gasout="gas.out"):
        # return E_solv-E_gas
        E_solv = self.get_Esp_from_out(out=solvout)
        E_gas = self.get_Esp_from_out(out=gasout)
        if any(np.isinf(x) for x in [E_gas, E_solv]):
            Gsolv_corr = np.inf
        else:
            Gsolv_corr = E_solv-E_gas
        return Gsolv_corr



class CallTerachem:
    _instance = None

    def __new__(cls, wd, Setting=None):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, wd, command_d):
        self.wd = wd
        self.cmd = mydict(command_d)

    def wirte_input(self):
        epsilon = self.cmd.get_value("epsilon")
        coordfn = self.cmd.get_value("coordfn")
        c = self.cmd.get_value("c")
        m = self.cmd.get_value("m")
        functional = self.cmd.get_value("functional")
        basis = self.cmd.get_value("basis")
        if m != 1:
            functional = "u" + functional

        with open(join(self.wd, "tc.inp"), 'w') as f:
            f.write(f"coordinates             {coordfn}.xyz\n"
                    f"method         {functional}\n"
                    f"basis                   {basis}\n"
                    "dftgrid                 1\n"
                    "sphericalbasis          yes\n"
                    "guess                   sad\n"
                    "precision               mixed\n"
                    "genscrdata              true\n"
                    "threall                 1e-14\n"
                    "convthre                1e-4\n"
                    "xtol                    1e-8\n"
                    f"charge                  {c}\n"
                    f"spinmult                {m}\n"
                    )
        if epsilon is None:
            with open(join(self.wd, "tc.inp"), 'a+') as f:
                f.write("run                     energy\n")
        else:
            with open(join(self.wd, "tc.inp"), 'a+') as f:
                f.write("pcm			cosmo\n")
                f.write(f"epsilon			{epsilon}\n")
                f.write("run                     energy\n")

    def run(self):
        cmd = f"{os.environ['TC_BIN']} tc.inp"
        # write input
        self.wirte_input()
        # run
        with open(join(self.wd, "tc.out"), 'w') as stdout_file, open('tc_err.out', 'w') as stderr_file:
            p = subprocess.Popen(cmd, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        return rc


class CallsTDA:
    def __init__(self, wd, coordfn):
        self.wd = wd
        self.coordfn = coordfn

    def callml(self, sty):
        cmd = f"{os.environ['STDA_BIN']} -f {self.coordfn} -ml_feature -sty {str(sty)}"
        p = subprocess.Popen(cmd, stdout=open(join(self.wd, f'stda.out'), 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc


class CallMultiwfn:
    def __init__(self, wd, inputfn):
        self.wd = wd
        self.inputfn = inputfn

    def write_cartesian_molden(self):
        with open(join(self.wd, "cmd.in"), "w") as f:
            f.write("100\n"
                    "2\n"
                    "6\n"
                    "orca.molden\n"
                    "0\n"
                    "-10\n"
                    )
        cmd = f"{os.environ['MWFN_BIN']} {self.inputfn}"
        with open(join(self.wd, "cmd.in"), "r") as cmd_in_file, open(join(self.wd, "multiwfn.out"), 'w') as out_file:
            p = subprocess.Popen(cmd, stdin=cmd_in_file, stdout=out_file, stderr=subprocess.PIPE, shell=True)
            p.wait()
        rc = p.returncode
        return rc


class CallMolBar:
    def __init__(self, coordinate, ele_l):
        self.coordinate = coordinate
        self.ele_l = ele_l

    def get_mb(self):
        from molbar.barcode import get_molbar_from_coordinates
        mb = get_molbar_from_coordinates(self.coordinate, self.ele_l)
        return mb


class CallgXTB:
    @classmethod
    def default_settings(cls, xtbsetting,mpsetting):
        # use gfn2 setting to obtain .charg and .uhf
        """Dynamically generates default settings based on the `setting` class."""
        return {"method": xtbsetting.method,
                "coordfn": "coord",
                "optlevel": xtbsetting.optlevel,
                "charge": 0,
                "uhf": 0,
                "alpb": xtbsetting.alpb,
                "etemp": xtbsetting.struc_etemp,
                "nthreads": mpsetting.xtb_threads
                }

    def __init__(self, wd, command_d=None, custom_settings=None):
        self.wd = wd
        if command_d is not None:
            self.cmd = mydict(command_d)
        else:
            import aRST.script.initialization as storage
            command_d = self.default_settings(storage.setting.xtb_setting,storage.setting.mp_setting)
            if custom_settings is not None:
                command_d.update(custom_settings)

            self.cmd = mydict(command_d)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(self.cmd["nthreads"])

    def gxtbsubprocess(self, cmd_string):
        out_path = f'{self.wd}/gxtb.out'
        err_path = f'{self.wd}/gxtb_err.out'
        with open(out_path, 'w') as stdout_file, open(err_path, 'w') as stderr_file:
            p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
            if rc == 0 and os.path.exists(err_path) and os.path.getsize(err_path) == 0:
                os.remove(err_path)
        return rc

    def generate_gXTBcmd_string(self):
        cmd_string = os.environ['GXTB_BIN'] + " -c "
        coordfn = self.cmd.get_value("coordfn")
        if coordfn is None:
            raise Exception(f"coord in {self.wd} can't find right coordfn name")
        cmd_string += f"{coordfn}.xyz "
        return cmd_string

    def sp(self):
        # write .CHRG and .UHF file
        if self.cmd.get_value('charge') != 0:
            with open(join(self.wd, ".CHRG"), 'w') as file:
                file.write(str(self.cmd.get_value('charge')))
        if self.cmd.get_value('uhf') != 0:
            with open(join(self.wd, ".UHF"), 'w') as file:
                file.write(str(self.cmd.get_value('uhf')))

        # run gxtb sp
        cmd_string2 = self.generate_gXTBcmd_string()
        rc = self.gxtbsubprocess(cmd_string2)
        return rc

    def afo_molden(self):
        if self.cmd.get_value('charge') != 0:
            with open(join(self.wd, ".CHRG"), 'w') as file:
                file.write(str(self.cmd.get_value('charge')))
        if self.cmd.get_value('uhf') != 0:
            with open(join(self.wd, ".UHF"), 'w') as file:
                file.write(str(self.cmd.get_value('uhf')))

        # run gxtb sp
        cmd_string2 = self.generate_gXTBcmd_string() + " -molden "
        rc = self.gxtbsubprocess(cmd_string2)
        return rc

    def get_Esp_from_out(self, out="gxtb.out"):
        import numpy as np
        # unit in eh
        Esp = np.inf
        if exists(join(self.wd, "NOT_CONVERGED")):
            return Esp

        if exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                lines = file.readlines()

            for ln in reversed(lines):
                if "total" in ln and not "time" in ln:
                    Esp = float(ln.split()[-1])
                    break
        return Esp


class CallFFiTS:
    @classmethod
    def default_settings(cls, xtbsetting):
        """Dynamically generates default settings based on the `setting` class."""
        return {
            "reactantfn": "reactant",
            "productfn": "product",
            "charge": 0,
            "uhf": 0,
            "alpb": xtbsetting.alpb,
            "etemp": xtbsetting.struc_etemp,
        }

    def __init__(self, wd, command_d=None, custom_settings=None):
        self.wd = wd
        if command_d is not None:
            self.cmd = mydict(command_d)
        else:
            import aRST.script.initialization as storage
            command_d = self.default_settings(storage.setting.xtb_setting)
            if custom_settings is not None:
                command_d.update(custom_settings)

            self.cmd = mydict(command_d)

    def xtbsubprocess(self, cmd_string):
        out_path = f'{self.wd}/xtb.out'
        err_path = f'{self.wd}/xtb_err.out'
        with open(out_path, 'w') as stdout_file, open(err_path, 'w') as stderr_file:
            p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        return rc

    def write_input(self):
        with open(join(self.wd, "input.toml"), "w") as f:
            f.write(f"""
[system]
charge={self.cmd.get('charge')}
multiplicity={self.cmd.get('uhf')}          
[reactant.path]
hessian_filename="rhessian"

[product.path]
hessian_filename="phessian"            
            """)

    @staticmethod
    def remove_trailing_blank_lines(filepath):
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()

        # remove trailing newlines at end of file
        content = content.rstrip("\r\n")

        with open(filepath, "w", encoding="utf-8") as f:
            f.write(content)

    def do(self):
        # do R hess
        cmd_string = f"{os.environ['XTB_BIN']} {self.cmd.get('reactantfn')}.xyz --hess -c {self.cmd.get('charge')} --uhf {self.cmd.get('uhf')} --alpb {self.cmd.get('alpb')} --etemp {self.cmd.get('etemp')}"
        rc = self.xtbsubprocess(cmd_string)
        if rc == 0:
            os.rename(join(self.wd, "hessian"), join(self.wd, "rhessian"))
            # do P hess
            cmd_string = f"{os.environ['XTB_BIN']} {self.cmd.get('productfn')}.xyz --hess -c {self.cmd.get('charge')} --uhf {self.cmd.get('uhf')} --alpb {self.cmd.get('alpb')} --etemp {self.cmd.get('etemp')}"
            rc = self.xtbsubprocess(cmd_string)
            if rc == 0:
                os.rename(join(self.wd, "hessian"), join(self.wd, "phessian"))

                # write input
                self.write_input()

                # # end-of-file error
                # self.remove_trailing_blank_lines(join(self.wd,self.cmd.get('reactantfn')))
                # self.remove_trailing_blank_lines(join(self.wd, self.cmd.get('productfn')))

                # run ffits
                out_path = f'{self.wd}/ffits.out'
                err_path = f'{self.wd}/ffits_err.out'
                cmd_string = f"ffits {self.cmd.get('reactantfn')} {self.cmd.get('productfn')} --config input.toml"
                with open(out_path, 'w') as stdout_file, open(err_path, 'w') as stderr_file:
                    p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
                    p.wait()
                    rc = p.returncode
        return rc