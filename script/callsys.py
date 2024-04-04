import builtins
import os
import shutil
import subprocess
import sys
import time
from tempfile import mkdtemp
import re
from func_timeout import func_set_timeout


# def print(text):
#     builtins.print(text)
#     os.fsync(sys.stdout)


def remove_file(old_path, new_path):
    filelist = os.listdir(old_path)

    for file in filelist:
        src = os.path.join(old_path, file)
        dst = os.path.join(new_path, file)

        shutil.move(src, dst)


class CallORCA:
    def __init__(self,level,coordfn,charge, m, pcm):
        self.coordfn = coordfn
        self.charge = charge
        self.m = m
        self.pcm = pcm
        self.level = level

    @staticmethod
    def orcasubprocess(cmd):
        rc = 1
        try:
            with open('orca.out', 'w') as stdout_file, open('orca_err.out', 'w') as stderr_file:
                p = subprocess.Popen(cmd, stdout=stdout_file, stderr=stderr_file, shell=True)
                p.wait()
                rc = p.returncode

                if rc == 0:
                    with open('orca.out', 'r') as f:
                        lines = f.readlines()
                        for line in lines:
                            if 'FINAL SINGLE' in line:
                                with open('ENERGY', 'w') as f2:
                                    f2.write(str(re.findall(r'\-?\d+\.?\d*', line)[0]))
                else:
                    print(f"Error: orca execution failed with return code {rc} at {os.getcwd()}")

        except Exception as e:
            print(f"An error occurred: {e}")
        return rc

    def callcalculation(self, caltyp='sp'):
        if self.m == 1:
            caltyp2 = "RKS"
        else:
            caltyp2 = "UKS"
        if self.pcm is not None:
            solmode = f"CPCM({self.pcm})"
        else:
            solmode = ''

        input = f"! {self.level} {caltyp} {caltyp2} {solmode}\n" \
                f"%SCF guess hueckel end\n" \
                f"* xyzfile {self.charge} {self.m} {self.coordfn}.xyz "
        with open(f'orca.inp', 'w') as f:
            f.write(input)
        cmd = f'$ORCA orca.inp'
        rc = self.orcasubprocess(cmd)

        return rc


class CallxTB:
    def __init__(self, level, coordfn, charge, uhf, alpb, etemp):
        self.coordfn = coordfn
        self.charge = charge
        self.uhf = uhf
        self.alpb = alpb
        self.etemp = etemp
        self.level = level

    @staticmethod
    def xtbsubprocess(cmd):
        rc = 1
        try:
            with open('xtb.out', 'w') as stdout_file, open('xtb_err.out', 'w') as stderr_file:
                p = subprocess.Popen(cmd, stdout=stdout_file, stderr=stderr_file, shell=True)
                p.wait()
                rc = p.returncode

                if rc == 0:
                    with open('xtb.out', 'r') as f:
                        lines = f.readlines()
                        for line in lines:
                            if 'TOTAL ENERGY' in line:
                                with open('ENERGY', 'w') as f2:
                                    f2.write(str(re.findall(r'\-?\d+\.?\d*', line)[0]))
                else:
                    print(f"Error: xtb execution failed with return code {rc}")

        except Exception as e:
            print(f"An error occurred: {e}")
        return rc

    def callsp(self):
        cmd = f'xtb {self.coordfn}.xyz --{self.level} -c {self.charge} --etemp {self.etemp}'

        if self.uhf is not None and self.uhf!=0:
            cmd += f' --uhf {self.uhf}'

        if self.alpb is None:
            cmd += f' --alpb {self.alpb}'

        rc = self.xtbsubprocess(cmd)
        return rc

    def callopt(self, extracmd=None):
        cmd = f'xtb {self.coordfn}.xyz --{self.level} -c {self.charge} --etemp {self.etemp}'

        if self.uhf is not None and self.uhf!=0:
            cmd += f' --uhf {self.uhf}'

        if self.alpb is not None:
            cmd += f' --alpb {self.alpb}'
        if extracmd is not None:
            cmd+= f' --opt {extracmd}'
        else:
            cmd += f' --opt '

        rc = self.xtbsubprocess(cmd)
        return rc
    def callConstrainedOpt(self,setting_d,extracmd=None):
        atl = setting_d["distance"]
        with open('xtb.inp', 'w') as f:
            f.write(f'$constrain\n')

            for i in range(len(atl)):
                f.write(f'distance: {atl[i][0]},{atl[i][1]},{atl[i][3]}\n')
            f.write("$end")
        cmd = f'xtb {self.coordfn}.xyz --{self.level} --input xtb.inp' \
              f' -c {self.charge} --uhf {self.uhf} --etemp {self.etemp}'
        if extracmd is not None:
            cmd+= f' --opt {extracmd}'
        else:
            cmd += f' --opt '
        p = subprocess.Popen(cmd, stdout=open(f'xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc



    def callml(self, xtbwd,mllevel=None,mletemp=None):
        cmd = f"{xtbwd} {self.coordfn}.xyz  -c {self.charge} --uhf {self.uhf} --ml_feature"
        if self.alpb is not None:
            cmd += f' --alpb {self.alpb}'

        if mllevel is not None:
            cmd+= f' --{mllevel}'
        if mletemp is not None:
            cmd += f' --etemp {mletemp}'

        p = subprocess.Popen(cmd, stdout=open(f'xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode

        return rc

    def calltrv(self):
        if self.alpb is None:
            command1 = f'xtb {self.coordfn}.xyz  --ohess --uhf {self.uhf} -c {self.charge}'
        else:
            command1 = f'xtb {self.coordfn}.xyz  --ohess --uhf {self.uhf} -c {self.charge} --alpb {self.alpb}'

        command2 = f'TRVthermo --io xtb --rotthr 50 --fxyz {self.coordfn}.xyz --fvib freq.out --temp 298.15 > thermo.out'
        p1 = subprocess.Popen(command1, stdout=open(f'freq.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p1.wait()
        p2 = subprocess.Popen(command2, shell=True)
        p2.wait()

        rc = p2.returncode
        return rc

    def callsolv(self):
        command1 = f'xtb {self.coordfn}.xyz --gfn 2 --uhf {self.uhf} -c {self.charge} bar1M '

        command2 = f'xtb {self.coordfn}.xyz --gfn 2 --alpb {self.alpb} --grad -- uhf {self.uhf} -c {self.charge} '
        p1 = subprocess.Popen(command1, stdout=open(f'gas.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p1.wait()
        rc1 = p1.returncode
        with open(f'gas.out', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'TOTAL ENERGY' in line:
                    egas = float(re.findall(r'\-?\d+\.?\d*', line)[0])
        p2 = subprocess.Popen(command2, stdout=open(f'solv.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p2.wait()
        with open(f'solv.out', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'TOTAL ENERGY' in line:
                    esolv = float(re.findall(r'\-?\d+\.?\d*', line)[0])
        rc2 = p2.returncode
        if rc1 == 0 and rc2 == 0:
            with open('GSOLV', 'w') as f2:
                f2.write(str(esolv - egas))

        return rc2

    def call_scan(self, setting_d):
        constrain_strength = setting_d["constrain"]["k"]
        constrain_distance = setting_d["constrain"]["distance"]
        mode = setting_d["scan"]["mode"]
        atl = setting_d["scan"]["distance"]

        step = setting_d["scan"]["step"]

        with open('scan.inp', 'w') as f:
            f.write(f'$constrain\n'
                    f'force constant={constrain_strength}\n')

            if constrain_distance is not None:
                f.write(f'distance:{constrain_distance[0]}, {constrain_distance[1]}, {constrain_distance[2]}\n')
            f.write(f'$scan\n')
            if mode is not None:
                f.write('mode=concerted\n')
            for i in range(len(atl)):
                f.write(f'distance: {atl[i][0]}, {atl[i][1]}, {atl[i][2]}; {atl[i][2]}, {atl[i][3]}, {step}\n')
            f.write("$end")
        cmd = f'xtb {self.coordfn}.xyz --{self.level} --input scan.inp --opt ' \
              f' -c {self.charge} --uhf {self.uhf} --etemp {self.etemp}'

        if self.alpb is not None:
            cmd += f' --alpb {self.alpb}'

        p = subprocess.Popen(cmd, stdout=open(f'xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        return rc
    def get_moldeninput(self):
        cmd = f'xtb {self.coordfn}.xyz --{self.level} -c {self.charge} --uhf {self.uhf} --etemp {self.etemp} --molden'

        if self.alpb is not None:
            cmd += f' --alpb {self.alpb}'


        rc = self.xtbsubprocess(cmd)
        return rc


class CallTerachem:
    def __init__(self,coordfn,method, charge, m):
        self.coordfn = coordfn
        self.method = method
        self.charge = charge
        self.m = m

    @staticmethod
    def tcsubprocess(cmd):
        rc = 1
        try:
            with open('tc.out', 'w') as stdout_file, open('tc_err.out', 'w') as stderr_file:
                p = subprocess.Popen(cmd, stdout=stdout_file, stderr=stderr_file, shell=True)
                p.wait()
                rc = p.returncode

        except Exception as e:
            print(f"An error occurred: {e}")
        return rc
    def get_moldeninput(self):
        with open(f'tc.inp', 'w') as f:
            f.write(f"coordinates             {self.coordfn}.xyz\n" 
                f"method         {self.method}\n" 
                "basis                   def2-sv_p\n" 
                "dftgrid                 1\n"
                "sphericalbasis          yes\n"
                "guess                   sad\n"
                "precision               mixed\n"
                "genscrdata              true\n"
                "threall                 1e-14\n"
                "convthre                1e-4\n"
                "xtol                    1e-6\n"    
                "pcm			cosmo\n"
                "epsilon			9.08\n"
                f"charge                  {self.charge}\n"
                f"spinmult                {self.m}\n"
                "run                     energy\n"
            )


        cmd = f'terachem tc.inp'
        rc = self.tcsubprocess(cmd)
        if rc ==0 and os.path.isdir(f"scr.{self.coordfn}"):
            cwd = os.getcwd()
            return os.path.join(cwd,f"scr.{self.coordfn}",f"{self.coordfn}.molden")
        else:
            return None


class CallsTDA:
    def __init__(self,coordfn):
        self.coordfn= coordfn

    @staticmethod
    def stdasubprocess(cmd):
        rc = 1
        try:
            with open('stda.out', 'w') as stdout_file, open('stda_err.out', 'w') as stderr_file:
                p = subprocess.Popen(cmd, stdout=stdout_file, stderr=stderr_file, shell=True)
                p.wait()
                rc = p.returncode

        except Exception as e:
            print(f"An error occurred: {e}")
        return rc
    def callml(self, stdawd,sty):
        cmd = f"{stdawd} -f {self.coordfn} -ml_feature -sty {str(sty)}"
        print(cmd)
        p = subprocess.Popen(cmd, stdout=open(f'stda.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode

        return rc









