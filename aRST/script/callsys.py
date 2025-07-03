import traceback
from os.path import join, exists
import os
import shutil
import subprocess
import re
from aRST.script.toolbox import mydict
import pandas as pd


def remove_file(old_path, new_path):
    filelist = os.listdir(old_path)

    for file in filelist:
        src = join(old_path, file)
        dst = join(new_path, file)
        shutil.move(src, dst)


class CallORCA:
    @classmethod
    def default_settings(cls, setting):
        """Dynamically generates default settings based on the `setting` class."""
        return {
            "functional": setting.get_value("orca_setting.functional", "r2scan-3c"),
            "basis": setting.get_value("orca_setting.basis", " "),
            "charge": 0,
            "m": 1,
            "fn": "coord",
            "caltyp": "SP",
            "caltyp2": "RKS",  # This can be updated dynamically if needed
            "solmode": setting.get_value("orca_setting.solmode", ""),
            "ncpu": setting.get_value("orca_setting.ncpu", 1),
        }

    def __init__(self, wd, command_d=None, custom_settings=None):
        self.wd = wd
        if command_d is not None:
            self.cmd = mydict(command_d)
        else:
            import aRST.script.initialization as storage
            command_d = self.default_settings(storage.setting)
            if custom_settings is not None:
                command_d.update(custom_settings)

            command_d["caltyp2"] = "RKS" if command_d["m"] == 1 else "UKS"
            self.cmd = mydict(command_d)

    def orcasubprocess(self, cmd_string, writeenergy=True, removegbw=True):
        with open(f'{self.wd}/orca.out', 'w') as stdout_file, open(f'{self.wd}/orca_err.out', 'w') as stderr_file:
            p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        if rc == 0 and writeenergy:
            with open(f'{self.wd}/orca.out', 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'FINAL SINGLE' in line:
                        with open(f'{self.wd}/ENERGY', 'w') as f2:
                            f2.write(str(re.findall(r'\-?\d+\.?\d*', line)[0]))
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

    def generate_cmd_string(self):
        command_string = "! "
        command_string += f"{self.cmd.get_value('functional', '')} "
        if self.cmd.get_value('basis', ''):
            command_string += f"{self.cmd.get_value('basis', '')} "
        command_string += f"{self.cmd.get_value('caltyp', '')} "
        command_string += f"{self.cmd.get_value('caltyp2', '')} "
        if self.cmd.get_value('solmode', ''):
            command_string += f"CPCM({self.cmd.get_value('solmode', '')}) \n"
        else:
            command_string += f"\n"

        ncpu = self.cmd.get_value('ncpu', 1)
        if isinstance(ncpu, int) and ncpu > 1:
            if self.cmd.get_value('ncpu', 1):
                command_string += f"%PAL NPROCS {int(ncpu)} END \n"
        c = self.cmd.get_value('charge')
        m = self.cmd.get_value('m')
        fn = self.cmd.get_value('fn')
        if all(v is not None for v in [fn, c, m]):
            command_string += f"* xyzfile {c} {m} {fn}.xyz "
        else:
            raise Exception(f"{fn}.xyz in {self.wd} has wrong c,m,fn parameter with {[fn, c, m]}")
        return command_string

    def call(self, removegbw=True):
        command_string = self.generate_cmd_string()
        with open(f'{self.wd}/orca.inp', 'w') as f:
            f.write(command_string)

        cmd = f"{os.environ['ORCA_BIN']} orca.inp"
        rc = self.orcasubprocess(cmd, removegbw=removegbw)
        return rc

    def get_geomres(self):
        if exists(join(self.wd, 'orca.xyz')):
            coord = pd.read_table(join(self.wd, 'orca.xyz'), header=None, skiprows=2, sep='\s+')
        else:
            coord = pd.read_table(join(self.wd, 'coord.xyz'), header=None, skiprows=2, sep='\s+')
        return coord

    def get_energy(self, out="orca.out"):
        Esp = None
        if exists(join(self.wd, "ENERGY")):
            with open(join(self.wd, "ENERGY"), 'r') as file:
                Esp = float(file.readlines()[0])
        elif exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                lines = file.readlines()
            for ln in lines:
                if "FINAL SINGLE POINT ENERGY" in ln:
                    Esp = float(ln.split()[-1])
                    break
        return Esp

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


class CallCrest:
    @classmethod
    def default_settings(cls, setting):
        """Dynamically generates default settings based on the `setting` class."""
        return {"method": setting.get_value("xtb_setting.method", "gfn2"),
                "coordfn": "coord",
                "optlevel": setting.get_value("xtb_setting.optlevel", ""),
                "charge": 0,
                "uhf": 0,
                "alpb": setting.get_value("xtb_setting.alpb", "h2o"),
                "etemp": setting.get_value("xtb_setting.struc_etemp", "4000"),
                }

    def __init__(self, wd, command_d=None, custom_settings=None):
        self.wd = wd
        if command_d is not None:
            self.cmd = mydict(command_d)
        else:
            import aRST.script.initialization as storage
            command_d = self.default_settings(storage.setting)
            if custom_settings is not None:
                command_d.update(custom_settings)

            self.cmd = mydict(command_d)

    def crestsubprocess(self, cmd_string):
        with open(f'{self.wd}/crest.out', 'w') as stdout_file, open(f'{self.wd}/crest_err.out', 'w') as stderr_file:
            p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        return rc

    def generate_cmd_string(self):
        cmd_string = os.environ['CREST_BIN']+" "
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

    def read_update_coord(self,job):
        outwd = f'{self.wd}/crest.out'
        with open(outwd, "r") as file:
            lines = file.readlines()
        trigger = False
        startline, end = 0, 0
        for line in lines:
            if "structure    ΔE(kcal/mol)   Etot(Eh)" in line:
                startline = lines.index(line)
                trigger = True
            if trigger and "Wall Time Summary" in line:
                end = lines.index(line) - 3
                trigger = False
                break
        deltaE_l = []
        for line in lines[startline+1:end]:
            nstruc,deltaE,Etot = line.split()
            deltaE_l.append(float(deltaE))
        if min(deltaE_l)>0:
            coord = None
        else:
            nstruc = deltaE_l.index(min(deltaE_l))
            ensemblewd = f'{self.wd}/deprotonated.xyz' if job=="depotonation" else f'{self.wd}/protonated.xyz'
            coord_l = []
            with open(ensemblewd, 'r') as f:

                while True:
                    # Read frame header
                    line = f.readline()
                    if not line:  # End of file
                        break
                    num_atoms = int(line.strip())
                    # Read comment line (contains energy)
                    comment_line = f.readline().strip()

                    # Read atomic coordinates
                    data = []
                    for _ in range(num_atoms):
                        parts = f.readline().split()
                        data.append([parts[0], float(parts[1]), float(parts[2]), float(parts[3])])

                    coord = pd.DataFrame(data)
                    coord_l.append(coord)
            coord = coord_l[nstruc]
        return coord




class CallxTB:
    @classmethod
    def default_settings(cls, setting):
        """Dynamically generates default settings based on the `setting` class."""
        return {"method": setting.get_value("xtb_setting.method", "gfn2"),
                "coordfn": "coord",
                "optlevel": setting.get_value("xtb_setting.optlevel", ""),
                "charge": 0,
                "uhf": 0,
                "alpb": setting.get_value("xtb_setting.alpb", "h2o"),
                "etemp": setting.get_value("xtb_setting.struc_etemp", "4000"),
                }

    def __init__(self, wd, command_d=None, custom_settings=None):
        self.wd = wd
        if command_d is not None:
            self.cmd = mydict(command_d)
        else:
            import aRST.script.initialization as storage
            command_d = self.default_settings(storage.setting)
            if custom_settings is not None:
                command_d.update(custom_settings)

            self.cmd = mydict(command_d)

    def xtbsubprocess(self, cmd_string, writeenergy=False):
        with open(f'{self.wd}/xtb.out', 'w') as stdout_file, open(f'{self.wd}/xtb_err.out', 'w') as stderr_file:
            p = subprocess.Popen(cmd_string, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        if rc == 0 and writeenergy:
            with open(f'{self.wd}/xtb.out', 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'TOTAL ENERGY' in line:
                        with open(f'{self.wd}/ENERGY', 'w') as f2:
                            f2.write(str(re.findall(r'\-?\d+\.?\d*', line)[0]))
        return rc

    def generate_cmd_string(self, opt=False, xtbmlversion=False,bar1M=False):
        if xtbmlversion:
            cmd_string = os.environ['XTBML_BIN'] + " "
        else:
            cmd_string = os.environ['XTB_BIN'] + " "
        coordfn = self.cmd.get_value("coordfn")
        if coordfn is None:
            raise Exception(f"coord in {self.wd} can't find right coordfn name")
        cmd_string += f"{coordfn}.xyz "
        if opt:
            optlevel = self.cmd.get_value("optlevel")
            cmd_string += f"--opt {optlevel} " if optlevel else "--opt "

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
        if bar1M and self.cmd.get_value('alpb'):
            cmd_string += "bar1M "
        cmd_string += f"{self.cmd.get_value('extracmd')} " \
            if self.cmd.get_value('extracmd') else ""
        return cmd_string

    def sp(self, writeenergy=True):
        cmd_string = self.generate_cmd_string(opt=False)
        rc = self.xtbsubprocess(cmd_string, writeenergy)
        return rc

    def opt(self, writeenergy=True):
        cmd_string = self.generate_cmd_string(opt=True)
        rc = self.xtbsubprocess(cmd_string, writeenergy)
        # if rc!=0:
        #     rc = self.secondchance_opt()
        return rc

    def afo(self):
        cmd_string = self.generate_cmd_string(opt=False, xtbmlversion=True)
        cmd_string += "--ml_feature "
        p = subprocess.Popen(cmd_string, stdout=open(f'{self.wd}/xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc

    def trv(self):
        cmd1 = self.generate_cmd_string()
        cmd1 += "--ohess "
        p1 = subprocess.Popen(cmd1, stdout=open(f'{self.wd}/freq.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p1.wait()
        rc = p1.returncode
        if rc == 0:
            cmd2 = f'{os.environ["TRVT_BIN"]} --io xtb --rotthr 50 '
            cmd2 += f'--fxyz xtbopt.xyz' if os.path.exists(
                join(self.wd, "xtbopt.xyz")) else f'--fxyz {self.cmd.get_value("coordfn")}.xyz '
            cmd2 += f" --temp {self.cmd.get_value('etemp')} " \
                if self.cmd.get_value('etemp') else ""
            cmd2 += f' --fvib freq.out > thermo.out'
            p2 = subprocess.Popen(cmd2, shell=True)
            p2.wait()
            rc2 = p2.returncode
            # if rc2 != 0:
            # raise Exception("TRVthermo calculation failed in",self.wd)
        else:
            rc2 = 1
            # raise Exception("TRV ohess calculation failed in ",self.wd)
        return rc2

    def solv(self):
        cmd2 = self.generate_cmd_string(bar1M=True)
        p2 = subprocess.Popen(cmd2, stdout=open(f'{self.wd}/solv.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p2.wait()
        with open(f'solv.out', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'TOTAL ENERGY' in line:
                    esolv = float(re.findall(r'\-?\d+\.?\d*', line)[0])
        rc2 = p2.returncode

        oldvalue = self.cmd.get_value("alpb")
        self.cmd.set_value("alpb", False)
        cmd1 = self.generate_cmd_string()
        self.cmd.set_value("alpb", oldvalue)

        p1 = subprocess.Popen(cmd1, stdout=open(f'{self.wd}/gas.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p1.wait()
        rc1 = p1.returncode
        with open(f'{self.wd}/gas.out', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'TOTAL ENERGY' in line:
                    egas = float(re.findall(r'\-?\d+\.?\d*', line)[0])
        if rc1 == 0 and rc2 == 0:
            with open(f'{self.wd}/GSOLV', 'w') as f2:
                f2.write(str(esolv - egas))

        return rc2

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

        cmd = self.generate_cmd_string(opt=True)
        cmd += "--input scan.inp "

        p = subprocess.Popen(cmd, stdout=open(f'{self.wd}/xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        if rc != 0 and (os.path.exists(f'{self.wd}/xtbscan.log') and os.path.exists(f'{self.wd}/xtbopt.xyz')):
            rc = 0
        return rc

    def constrained_opt(self, newbond):
        with open(f'{self.wd}/scan.inp', 'w') as f:
            f.write(f'$constrain\n'
                    f'force constant={self.cmd.get_value("k","1")}\n')
            for item in newbond:
                at1, at2 = item["atoms"]
                value = item["value"]
                f.write('distance: %d,%d,%.2f\n'
                        % (at1, at2, value))
            f.write("$end")
        cmd = self.generate_cmd_string(opt=True)
        cmd += "--input scan.inp "
        p = subprocess.Popen(cmd, stdout=open(f'{self.wd}/xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        if rc != 0 and os.path.exists(f'{self.wd}/xtbopt.xyz'):
            rc = 0
        return rc

    def secondchance_opt(self):
        cmd_string = "xtb "
        coordfn = self.cmd.get_value("coordfn")
        cmd_string += f"{coordfn}.xyz "
        cmd_string += "--opt "

        cmd_string += f"-c {self.cmd.get_value('charge')} " \
            if self.cmd.get_value('charge') else ""
        cmd_string += f"--uhf {self.cmd.get_value('uhf')} " \
            if self.cmd.get_value('uhf') else ""

        rc = self.xtbsubprocess(cmd_string, True)
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
        cmd_string = f"xtb {self.cmd.get_value('coordfn')}.xyz --input md.inp --md "
        cmd_string += f"-c {self.cmd.get_value('charge')} " \
            if self.cmd.get_value('charge') else ""
        cmd_string += f"--uhf {self.cmd.get_value('uhf')} " \
            if self.cmd.get_value('uhf') else ""
        cmd_string += f"--alpb {self.cmd.get_value('alpb')} " \
            if self.cmd.get_value('alpb') else ""
        p = subprocess.Popen(cmd_string, stdout=open(f'{self.wd}/xtbmd.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc

    def get_Emin_geom_from_mdtrj(self, nat):
        import numpy as np

        # cut geomblock in md.trj
        coord_l = []
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

                    # Read atomic coordinates
                    data = []
                    for _ in range(num_atoms):
                        parts = f.readline().split()
                        data.append([parts[0], float(parts[1]), float(parts[2]), float(parts[3])])

                    coord = pd.DataFrame(data)
                    coord_l.append(coord)
                    energies.append(energy)
            return coord_l, energies

        def filter_outliers(coord_l, energies, z_threshold=2.5):
            from scipy.stats import zscore
            z_scores = np.abs(zscore(energies))
            non_outliers = np.where(z_scores < z_threshold)[0]

            clean_coord_l = [coord_l[i] for i in non_outliers]
            clean_energies = [energies[i] for i in non_outliers]

            return clean_coord_l, clean_energies

        frames, energies = parse_xtb_trj(trjwd, nat)
        clean_coord_l, clean_energies = filter_outliers(coord_l, energies)
        lowest_idx = np.argmin(clean_energies)
        coord_emin = clean_coord_l[lowest_idx]
        emin = clean_energies[lowest_idx]
        return coord_emin, emin

    def get_moldeninput(self):
        cmd = self.generate_cmd_string()
        cmd += "--molden "
        rc = self.xtbsubprocess(cmd)
        return rc

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
        Esp = None
        if exists(join(self.wd, "ENERGY")):
            with open(join(self.wd, "ENERGY"), 'r') as file:
                Esp = float(file.readlines()[0])

        elif exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                lines = file.readlines()

            for ln in lines:
                if "total energy" in ln:
                    Esp = float(ln.split()[4])
                    break
        return Esp

    def get_Etrv(self, out="GTRV"):
        # unit in eh
        Etrv = None
        if exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                Etrv = float(file.readlines()[0])
        return Etrv

    def get_Esolv(self, out="GSOLV"):
        # unit in eh
        Esolv = None
        if exists(join(self.wd, out)):
            with open(join(self.wd, out), 'r') as file:
                Esolv = float(file.readlines()[0])
        return Esolv

    def get_topo(self, out="xtbtopo.mol"):
        # constrain_l, ele required in mb; at index starts with 1
        constrain_l, ele = None, None
        coord = None
        for file in os.listdir(self.wd):
            if file.endswith("xyz"):
                from aRST.script.toolbox import read_coord
                coord = read_coord(join(self.wd, file))
                break
        if coord is not None:
            if exists(join(self.wd, out)):
                constrain_l, ele = [], []
                from aRST.script.geom import Coord
                # V2000
                with open(join(self.wd, out)) as f:
                    lines = f.readlines()
                ln_parts = lines[3].split()
                nat, nbond = int(ln_parts[0]), int(ln_parts[1])

                for i in range(nat):
                    ele.append(lines[4 + i].split()[3])
                for i in range(nbond):
                    ln_parts = lines[4 + nat + i].split()
                    at1, at2 = int(ln_parts[0]), int(ln_parts[1])
                    bond_length = Coord(coord).get_ats_distance(at1 - 1, at2 - 1)
                    constrain_l.append({"atoms": [at1, at2], "value": bond_length})
            else:
                # single at
                constrain_l, ele = [], [str(coord.iloc[0, 0])]
        return constrain_l, ele


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
    def __init__(self, coord):
        self.coord = coord

    def get_mb(self):
        from molbar.barcode import get_molbar_from_coordinates
        coord_l = [self.coord.loc[i, 1:3].tolist() for i in self.coord.index]
        ele_l = [self.coord.loc[i, 0] for i in self.coord.index]
        mb = get_molbar_from_coordinates(coord_l, ele_l)
        return mb
