import os
import shutil
import subprocess
import time
from tempfile import mkdtemp
import re


def remove_file(old_path, new_path):
    filelist = os.listdir(old_path)

    for file in filelist:
        src = os.path.join(old_path, file)
        dst = os.path.join(new_path, file)

        shutil.move(src, dst)


def call_orca_dft(wd, method, charg, m, nat, opt=True, alpb=None, mode=0):
    os.chdir(wd)
    if nat > 1:
        caltyp = 'sp OPT'
    else:
        caltyp = 'sp'
    if not opt:  # Manually disable optimization
        caltyp = 'sp'
    solmode = ""
    if charg:
        if alpb is None:
            print(f"please use CPCM mode for {wd}")
        else:
            if alpb == "h2o" or "H2O":
                alpb = "water"
            solmode = f"CPCM({alpb})"
    if m == 1:
        caltyp2 = "RKS"
    else:
        caltyp2 = "UKS"
    cmd = f"! {method} {caltyp} {caltyp2} {solmode}\n" \
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
        call_subprocess(cmd,'orca')
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
        call_subprocess(cmd,'o')
        return


def call_xtb(wd, charg, nat, alpb=None, mode=0):
    os.chdir(wd)
    if nat > 1:
        caltyp = '--opt'
    else:
        caltyp = ' '
    solmode = ''
    if charg:
        if alpb is None:
            print(f"please use alpb mode for {wd}")
        else:
            solmode = f'--alpb {alpb}'
    command = f'xtb coord.xyz {caltyp} -c {charg} {solmode}' \
              # f'&& tac xtb.out | grep -a \'TOTAL ENERGY\' |awk \'{{print $4}}\' > ENERGY &\n'

    if mode == 1:  # running in /tmp and return when cal has finished
        tmpd = mkdtemp()
        remove_file(wd, tmpd)
        os.chdir(tmpd)
        call_subprocess(command,'xtb')
        remove_file(tmpd, wd)
        os.chdir(wd)
        if os.path.isdir(tmpd):
            shutil.rmtree(tmpd)
        return
    elif mode == 2:  # running in master node but submit jobs using slurm
        root_path = os.path.dirname(__file__)
        shutil.copy(f'{root_path}/submit/submit_xtbopt.sh', wd)
        os.chdir(wd)
        cmd = f'sbatch submit_xtbopt.sh'
        subprocess.run(cmd, shell=True)
        return
    else:
        call_subprocess(command,'xtb')
        return


def call_xtb_scan(wd, charg, alpb=None, mode=0, **kw):
    os.chdir(wd)
    concerted_trigger = kw.get('concerted_trigger', False)
    etmp = kw.get('etmp', str(300))
    kw = kw.get('kw')
    at1 = kw['at1']
    at2 = kw['at2']
    d1 = kw['d1']
    d2 = kw['d2']

    if concerted_trigger:
        at3 = kw['at3']
        at4 = kw['at4']
        d_34_start = kw['d_34_start']
        d_34_end = kw['d_34_end']
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant=1\n'
                     f'$scan\n'
                     'mode=concerted\n'
                     f'distance: {at3},{at4},{d_34_start};{d_34_start},{d_34_end},30\n'
                     f'distance: {at1},{at2},{d1};{d1},{d2},30\n$end')
    else:
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant=1\n'
                     f'$scan\n'
                     f'distance: {at1},{at2},{d1};{d1},{d2},30\n$end')
    solmode = ""
    if charg:
        if alpb is None:
            print(f"please use alpb mode for {wd}")
        else:
            solmode = f"--alpb {alpb}"
    command = f'xtb coord.xyz --opt  -c {charg}  {solmode} --etemp {etmp} --input scan.inp'

    if mode == 1:
        tmpd = mkdtemp()
        remove_file(wd, tmpd)
        os.chdir(tmpd)
        call_subprocess(command,'xtb')
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
        call_subprocess(command, 'xtb')
        return


def call_subprocess(cmd,typ):
    keyw = None
    if typ == 'orca':
        keyw = 'FINAL SINGLE'
    if typ == 'xtb':
        keyw = 'TOTAL ENERGY'
    if keyw:
        p = subprocess.Popen(cmd, stdout=open(f'{typ}.out', 'w'), stderr=subprocess.PIPE,shell=True)
        p.wait()
        with open(f'{typ}.out', 'r') as f:
            lines = f.readlines()
            line = re.findall(r'\-?\d+\.?\d*', [k for k in lines if keyw in k][0])
            e = float('inf')
            if line:
                try:
                    e = float(line[0])
                except:
                    pass
        with open('ENERGY', 'w') as f:
            f.write(str(e))
    return
