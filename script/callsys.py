import builtins
import os
import shutil
import subprocess
import sys
import time
from tempfile import mkdtemp
import re
from func_timeout import func_set_timeout


def print(text):
    builtins.print(text)
    os.fsync(sys.stdout)


def remove_file(old_path, new_path):
    filelist = os.listdir(old_path)

    for file in filelist:
        src = os.path.join(old_path, file)
        dst = os.path.join(new_path, file)

        shutil.move(src, dst)


def call_orca_dft(wd, charg, m, caltyp='sp', **kwargs):
    os.chdir(wd)
    pcm = kwargs.get("pcm", None)
    method = kwargs.get("dft_level", "r2scan-3c")
    mode = kwargs.get("mode", 0)

    if pcm is not None:
        solmode = f"CPCM({pcm})"
    else:
        solmode = ''

    if m == 1:
        caltyp2 = "RKS"
    else:
        caltyp2 = "UKS"
    # f"pal nprocs 16 end\n"\
    # f"%SCF MAXITER 200 end\n" \
    if caltyp == 'sp':
        # sp calculation only used for initial evaluation for reaction path
        cmd = f"! {method} {caltyp} {caltyp2} {solmode}\n" \
              f"%SCF guess hueckel end\n" \
              f"* xyzfile {charg} {m} coord.xyz "
    else:
        cmd = f"! {method} {caltyp} {caltyp2} {solmode}\n" \
              f"%SCF guess hueckel end\n" \
              f"* xyzfile {charg} {m} coord.xyz "
    with open(f'{wd}/orca.inp', 'w') as f:
        f.write(cmd)
    cmd = f'$ORCA orca.inp'
    # 'orca.out'
    #  && tac orca.out | grep -m 1 \'FINAL SINGLE\' |awk \'{{print $5}}\' > ENERGY
    if mode == 1:  # running in debug node and return when cal has finished
        tmpd = mkdtemp()
        remove_file(wd, tmpd)
        os.chdir(tmpd)
        rc = call_subprocess(cmd, 'orca')
        if ".tmp" in os.listdir(tmpd):
            cmd = f"! {method} sp OPT {caltyp2} {solmode}\n" \
                  f"" \
                  f"* xyzfile {charg} {m} coord.xyz "
            rc = call_subprocess(cmd, 'orca')
        if rc != 0:
            print(f"{wd} calling orca wrong")

        remove_file(tmpd, wd)
        os.chdir(wd)
        if os.path.isdir(tmpd):
            shutil.rmtree(tmpd)
        return

    elif mode == 2:  # running in master node but submit jobs using slurm
        root_path = os.path.dirname(__file__)
        shutil.copy(f'{root_path}/submit/submit_orca.sh', wd)
        os.chdir(wd)
        cmd = f'sbatch submit_orca.sh'
        subprocess.run(cmd, shell=True)
        return
    else:  # default
        rc = call_subprocess(cmd, 'orca')
        if rc != 0:
            print(f"{wd} calling orca wrong")
        return rc


def call_xtb(wd, charg, caltyp="", uhf=0, **kwargs):
    os.chdir(wd)
    alpb = kwargs.get("alpb", None)
    mode = kwargs.get("mode", 0)
    etemp = kwargs.get("etemp", str(300))
    if alpb is None:
        solmode = ""
    else:
        solmode = f'--alpb {alpb}'
    command = f'xtb coord.xyz {caltyp} -c {charg} --uhf {uhf} {solmode} --etemp {etemp}'

    # f'&& tac xtb.out | grep -a \'TOTAL ENERGY\' |awk \'{{print $4}}\' > ENERGY &\n'
    # if uhf is not None:
    #     command = f'xtb coord.xyz --opt --uhf {uhf} {solmode}'
    if mode == 1:  # running in /tmp and return when cal has finished
        tmpd = mkdtemp()
        remove_file(wd, tmpd)
        os.chdir(tmpd)
        rc = call_subprocess(command, 'xtb')
        if rc != 0:
            print("calling xtb wrong")
        remove_file(tmpd, wd)
        os.chdir(wd)
        if os.path.isdir(tmpd):
            shutil.rmtree(tmpd)
        return rc
    elif mode == 2:  # running in master node but submit jobs using slurm
        root_path = os.path.dirname(__file__)
        shutil.copy(f'{root_path}/submit/submit_xtbopt.sh', wd)
        os.chdir(wd)
        cmd = f'sbatch submit_xtbopt.sh'
        subprocess.run(cmd, shell=True)
        return
    else:
        rc = call_subprocess(command, 'xtb')
        if rc != 0:
            print("calling xtb wrong")
        return rc


def call_xtb_scan(wd, charg, uhf, **kwargs):
    os.chdir(wd)
    at1 = kwargs.get('at1')
    at2 = kwargs.get('at2')
    d1 = kwargs.get('d1')
    d2 = kwargs.get('d2')
    alpb = kwargs.get("alpb", None)
    mode = kwargs.get("mode", 0)
    etemp = kwargs.get("etemp", str(300))
    concerted_trigger = kwargs.get("concerted_trigger", False)
    constrain_strength = kwargs.get("constrain_strength", str(1))
    at3 = kwargs.get('at3', None)
    at4 = kwargs.get('at4', None)
    d_34_start = kwargs.get('d_34_start', None)
    d_34_end = kwargs.get("d_34_end", None)

    if concerted_trigger:
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant={constrain_strength}\n'
                     f'$scan\n'
                     'mode=concerted\n'
                     f'distance: {at3},{at4},{d_34_start};{d_34_start},{d_34_end},30\n'
                     f'distance: {at1},{at2},{d1};{d1},{d2},30\n$end')
    else:
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant={constrain_strength}\n'
                     f'$scan\n'
                     f'distance: {at1},{at2},{d1};{d1},{d2},30\n$end')

    if alpb is None:
        solmode = ""
    else:
        solmode = f'--alpb {alpb}'
    command = f'xtb coord.xyz --opt  -c {charg} --uhf {uhf} {solmode} --etemp {etemp} --input scan.inp'

    if mode == 1:
        tmpd = mkdtemp()
        remove_file(wd, tmpd)
        os.chdir(tmpd)
        rc = call_subprocess(command, 'xtb')
        if rc != 0:
            print("calling scan wrong")
        remove_file(tmpd, wd)
        os.chdir(wd)
        if os.path.isdir(tmpd):
            shutil.rmtree(tmpd)
    elif mode == 2:
        root_path = os.path.dirname(__file__)
        shutil.copy(f'{root_path}/submit/submit_scan.sh', wd)
        with open('submit_scan.sh', 'r') as f:
            line = f.readlines()
        line[45] = f'{command}'
        with open('submit_scan.sh', 'w') as f:
            for i in line:
                f.write(i)
        subprocess.run("sbatch submit_scan.sh", shell=True)
        return
    else:
        rc = call_subprocess(command, 'xtb')
        if rc != 0:
            print("calling scan wrong")
        return


def call_xtb_ml(wd, charg, uhf=0, **kwargs):
    os.chdir(wd)
    alpb = kwargs.get("alpb", None)
    xtbwd = kwargs.get("xtb_wd")
    etemp = kwargs.get("etemp", str(300))
    if alpb is not None:
        cmd = f'{xtbwd} coord.xyz --ml_feature --etemp {etemp} -c {charg} --uhf {uhf} --alpb {alpb}'
    else:
        cmd = f'{xtbwd} coord.xyz --ml_feature --etemp {etemp} -c {charg} --uhf {uhf}'

    # run ml
    # subprocess.run(cmd, shell=True)
    p = subprocess.Popen(cmd, stdout=open(f'xtb.out', 'w'), stderr=subprocess.PIPE, shell=True)
    p.wait()
    return


def call_trv(wd, charg, uhf=0, **kwargs):
    os.chdir(wd)
    alpb = kwargs.get("alpb", None)
    mode = kwargs.get("mode", 0)
    if alpb is None:
        solmode = ""
    else:
        solmode = f'--alpb {alpb}'
    command1 = f'xtb coord.xyz --gfn 2 --ohess --uhf {uhf} -c {charg} {solmode}'
    command2 = f'TRVthermo --io xtb --rotthr 50 --fxyz coord.xyz --fvib freq.out --temp 298.15 > thermo.out'
    rc = 0
    if mode == 0:
        p1 = subprocess.Popen(command1, stdout=open(f'freq.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p1.wait()
        p2 = subprocess.Popen(command2, shell=True)
        p2.wait()

        rc = p2.returncode
    if rc != 0:
        print("calling trv wrong")
    return rc


def call_solv(wd, charg, uhf=0, **kwargs):
    os.chdir(wd)
    alpb = kwargs.get("alpb", None)
    mode = kwargs.get("mode", 0)
    esolv = 0
    egas = 0
    if alpb is None:
        solmode = ""
    else:
        solmode = f'--alpb {alpb}'
    command1 = f'xtb coord.xyz --gfn 2 --uhf {uhf} -c {charg} bar1M '

    command2 = f'xtb coord.xyz --gfn 2 {solmode} --grad -- uhf {uhf} -c {charg} '
    rc = 0
    if mode == 0:
        p1 = subprocess.Popen(command1, stdout=open(f'gas.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p1.wait()
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
        if rc ==0:
            with open('GSOLV', 'w') as f2:
                f2.write(str(esolv - egas))
        rc = p2.returncode
    if rc != 0:
        print("calling trv wrong")
    return rc


def call_subprocess(cmd, typ):
    # call&check program running status and get calculated energy
    keyw = None
    if typ == 'orca':
        keyw = 'FINAL SINGLE'
    if typ == 'xtb':
        keyw = 'TOTAL ENERGY'
    if keyw:
        p = subprocess.Popen(cmd, stdout=open(f'{typ}.out', 'w'), stderr=subprocess.PIPE, shell=True)
        p.wait()
        rc = p.returncode
        with open(f'{typ}.out', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if keyw in line:
                    with open('ENERGY', 'w') as f2:
                        f2.write(str(re.findall(r'\-?\d+\.?\d*', line)[0]))
    else:
        rc = 1
    return rc
