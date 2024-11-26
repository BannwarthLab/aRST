from os.path import join
import os
import shutil
import subprocess
import re
from aRST.script.dictionaryInspector import DictionaryInspector


def remove_file(old_path, new_path):
    filelist = os.listdir(old_path)

    for file in filelist:
        src = join(old_path, file)
        dst = join(new_path, file)

        shutil.move(src, dst)


class CallORCA:
    def __init__(self, wd, command_d):
        self.wd = wd
        self.command_d_inspertor = DictionaryInspector(command_d)

    def orcasubprocess(self, cmd_string, writeenergy=True,removegbw=True):
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
        extensions_to_remove = [ '.densities', '.opt']
        if removegbw:
            extensions_to_remove.append('.gbw')
        for filename in os.listdir(f'{self.wd}'):
            if any(filename.endswith(ext) for ext in extensions_to_remove):
                os.remove(filename)
        return rc

    def generate_cmd_string(self):
        command_string = "! "
        command_string += f"{self.command_d_inspertor.get_value('functional', '')} "
        if self.command_d_inspertor.get_value('basis', ''):
            command_string += f"{self.command_d_inspertor.get_value('basis', '')} "
        command_string += f"{self.command_d_inspertor.get_value('caltyp', '')} "
        command_string += f"{self.command_d_inspertor.get_value('caltyp2', '')} "
        if self.command_d_inspertor.get_value('solmode', ''):
            command_string += f"CPCM({self.command_d_inspertor.get_value('solmode', '')}) \n"
        else:command_string += f"\n"


        ncpu = self.command_d_inspertor.get_value('ncpu', 1)
        if isinstance(ncpu, int) and ncpu > 1:
            if self.command_d_inspertor.get_value('ncpu', 1):
                command_string += f"%PAL NPROCS {int(ncpu)} END \n"
        c = self.command_d_inspertor.get_value('charge')
        m = self.command_d_inspertor.get_value('m')
        fn = self.command_d_inspertor.get_value('fn')
        if all(v is not None for v in [fn, c, m]):
            command_string += f"* xyzfile {c} {m} {fn}.xyz "
        else:
            raise Exception(f"{fn}.xyz in {self.wd} has wrong c,m,fn parameter with {[fn, c, m]}")
        return command_string

    def call(self,removegbw=True):
        command_string = self.generate_cmd_string()
        with open(f'{self.wd}/orca.inp', 'w') as f:
            f.write(command_string)

        cmd = f'$ORCA orca.inp'
        rc = self.orcasubprocess(cmd,removegbw=removegbw)
        return rc


class CallxTB:
    def __init__(self, wd, command_d):
        self.wd = wd
        self.command_d_inspertor = DictionaryInspector(command_d)

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

    def generate_cmd_string(self, opt=False, xtbmlversion=False):
        if xtbmlversion:
            cmd_string = self.command_d_inspertor.get_value("xtbwd")+" "
        else:
            cmd_string = "xtb "
        coordfn = self.command_d_inspertor.get_value("coordfn")
        if coordfn is None:
            raise Exception(f"coord in {self.wd} can't find right coordfn name")
        cmd_string += f"{coordfn}.xyz "
        if opt:
            optlevel = self.command_d_inspertor.get_value("optlevel")
            cmd_string += f"--opt {optlevel} " if optlevel else "--opt "

        cmd_string += f"--{self.command_d_inspertor.get_value('method')} " \
            if self.command_d_inspertor.get_value('method') else ""
        cmd_string += f"-c {self.command_d_inspertor.get_value('charge')} " \
            if self.command_d_inspertor.get_value('charge') else ""
        cmd_string += f"--uhf {self.command_d_inspertor.get_value('uhf')} " \
            if self.command_d_inspertor.get_value('uhf') else ""
        cmd_string += f"--etemp {self.command_d_inspertor.get_value('etemp')} " \
            if self.command_d_inspertor.get_value('etemp') else ""
        cmd_string += f"--alpb {self.command_d_inspertor.get_value('alpb')} " \
            if self.command_d_inspertor.get_value('alpb') else ""
        cmd_string += f"{self.command_d_inspertor.get_value('extracmd')} " \
            if self.command_d_inspertor.get_value('extracmd') else ""
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

    def ml(self):
        cmd_string = self.generate_cmd_string(opt=False, xtbmlversion=True)
        cmd_string+= "--ml_feature "
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
            cmd2 = f'TRVthermo --io xtb --rotthr 50 '
            cmd2 += f'--fxyz xtbopt.xyz' if os.path.exists(
                join(self.wd,"xtbopt.xyz")) else f'--fxyz {self.command_d_inspertor.get_value("coordfn")}.xyz '
            cmd2 += f" --temp {self.command_d_inspertor.get_value('etemp')} " \
                if self.command_d_inspertor.get_value('etemp') else ""
            cmd2 += f' --fvib freq.out > thermo.out'
            p2 = subprocess.Popen(cmd2, shell=True)
            p2.wait()
            rc2 = p2.returncode
            # if rc2 != 0:
                # raise Exception("TRVthermo calculation failed in",self.wd)
        else:
            rc2=1
            # raise Exception("TRV ohess calculation failed in ",self.wd)
        return rc2

    def solv(self):
        cmd2 = self.generate_cmd_string()
        p2 = subprocess.Popen(cmd2, stdout=open(f'{self.wd}/solv.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p2.wait()
        with open(f'solv.out', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'TOTAL ENERGY' in line:
                    esolv = float(re.findall(r'\-?\d+\.?\d*', line)[0])
        rc2 = p2.returncode

        self.command_d_inspertor.set_value("alpb",False)
        cmd1 = self.generate_cmd_string()

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
                f.write('distance: %d,%d,%.2f\n' % (at1, at2, dd * 1.2))

            f.write('$scan\n')
            if setting_d["scan"]["mode"] == "concerted":
                f.write('mode=concerted\n')

            atl = setting_d["scan"]["distance"]
            for item in atl:
                at1,at2,d0,d_end = item
                f.write('distance: %d,%d,%.2f;%.2f,%.2f,%d\n'
                        % (at1,at2, d0, d0, d_end, setting_d["scan"]["step"]))
            f.write("$end")

        cmd = self.generate_cmd_string(opt=True)
        cmd+= "--input scan.inp "

        p = subprocess.Popen(cmd, stdout=open(f'{self.wd}/xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        if rc !=0 and (os.path.exists(f'{self.wd}/xtbscan.log') and os.path.exists(f'{self.wd}/xtbopt.xyz')):
            rc=0
        return rc

    def constrained_opt(self,newbond):
        with open(f'{self.wd}/scan.inp', 'w') as f:
            f.write(f'$constrain\n'
                    f'force constant={self.command_d_inspertor.get_value("k")}\n')
            for item in newbond:
                at1,at2 = item["atoms"]
                value = item["value"]
                f.write('distance: %d,%d,%.2f\n'
                        % (at1,at2,value))
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
        coordfn = self.command_d_inspertor.get_value("coordfn")
        cmd_string += f"{coordfn}.xyz "
        cmd_string += "--opt "

        cmd_string += f"-c {self.command_d_inspertor.get_value('charge')} " \
            if self.command_d_inspertor.get_value('charge') else ""
        cmd_string += f"--uhf {self.command_d_inspertor.get_value('uhf')} " \
            if self.command_d_inspertor.get_value('uhf') else ""

        rc = self.xtbsubprocess(cmd_string, True)
        return rc



    def get_moldeninput(self):
        cmd = self.generate_cmd_string()
        cmd += "--molden "
        rc = self.xtbsubprocess(cmd)
        return rc


class CallTerachem:
    def __init__(self, wd, command_d):
        self.wd = wd
        self.command_d_inspertor = DictionaryInspector(command_d)

    def wirte_input(self):
        epsilon = self.command_d_inspertor.get_value("epsilon")
        coordfn = self.command_d_inspertor.get_value("coordfn")
        c = self.command_d_inspertor.get_value("c")
        m = self.command_d_inspertor.get_value("m")
        functional = self.command_d_inspertor.get_value("functional")
        basis = self.command_d_inspertor.get_value("basis")
        if m !=1:
            functional = "u"+functional

        with open(join(self.wd,"tc.inp"), 'w') as f:
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
        cmd = "terachem tc.inp"
        # write input
        self.wirte_input()
        # run
        with open(join(self.wd,"tc.out"), 'w') as stdout_file, open('tc_err.out', 'w') as stderr_file:
            p = subprocess.Popen(cmd, stdout=stdout_file, stderr=stderr_file, shell=True)
            p.wait()
            rc = p.returncode
        return rc

class CallsTDA:
    def __init__(self,wd, coordfn):
        self.wd = wd
        self.coordfn = coordfn

    def callml(self, stdawd, sty):
        cmd = f"{stdawd} -f {self.coordfn} -ml_feature -sty {str(sty)}"
        p = subprocess.Popen(cmd, stdout=open(join(self.wd,f'stda.out'), 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc

class CallMultiwfn:
    def __init__(self, wd,inputfn):
        self.wd = wd
        self.inputfn = inputfn

    def write_cartesian_molden(self):
        with open(join(self.wd,"cmd.in"),"w") as f:
            f.write("100\n"
                    "2\n"
                    "6\n"
                    "orca.molden\n"
                    "0\n"
                    "-10\n"
                    )
        cmd = f"/home/guests/ychen/Multiwfn/Multiwfn_noGUI {self.inputfn}"
        with open(join(self.wd, "cmd.in"), "r") as cmd_in_file, open(join(self.wd,"multiwfn.out"), 'w') as out_file:
            p = subprocess.Popen(cmd, stdin=cmd_in_file,stdout=out_file, stderr=subprocess.PIPE, shell=True)
            p.wait()
        rc = p.returncode
        return rc



