import os
import re
import shutil
import string
import random
import numpy as np



def create_folder(wd):
    wd = os.path.abspath(wd)
    if not os.path.exists(wd) and not os.path.isdir(wd):
        try:
            os.makedirs(wd)
        except:
            pass
    os.chdir(wd)
    return wd

def rm_folderORfile(wd):
    wd = os.path.abspath(wd)
    if os.path.isdir(wd):
        try:
            shutil.rmtree(f'{wd}', ignore_errors=True)
            os.rmdir(f'{wd}')
        except:
            pass
    elif os.path.isfile(wd):
        try:
            os.remove(wd)
        except:
            pass


def get_neighborats(at,wbo,filetype="bo"):
    neighbor_l = []
    if filetype=="bo":
        boindex = 2
    else:
        boindex = 3
    for index in range(wbo.shape[0]):
        at1, at2 = wbo.iloc[index, :2]
        dgap = wbo.iloc[index, boindex]
        if filetype=="bo" and dgap<=0.3:
            continue
        elif filetype!="bo" and wbo.iloc[index, boindex]<=0.3:
            continue
        if int(at) ==int(at1):
            neighbor_l.append(int(at2))
        elif int(at) ==int(at2):
            neighbor_l.append(int(at1))
    return list(set(neighbor_l))

def add_fragmentat(at,wbo,check_l,atomlist = []):
    if wbo is None:
        # single atom
        atomlist.append(at)
        check_l[at-1] = True
        return check_l,atomlist
    else:
        trigger = True
        while trigger:
            neighbor_l = get_neighborats(at,wbo)
            if not neighbor_l:
                # single atom
                atomlist.append(at)
                check_l[at - 1] = True
                trigger=False
            else:
                if not all(check_l[i-1] for i in neighbor_l):
                    for i in neighbor_l:
                        if not check_l[i-1]:
                            check_l[i-1] = True
                            atomlist.append(i)
                            if i not in atomlist:
                                atomlist.append(i)
                            check_l,atomlist = add_fragmentat(i,wbo,check_l,atomlist)
                else:
                    trigger = False
    return check_l,atomlist

def get_fr_valenz_nel(coord):
    from aRST.script.element import valence_electrons
    nel = 0
    for at in range(coord.shape[0]):
        Z = valence_electrons[coord.iloc[at,0]]
        nel +=Z
    return nel

import math
def round_charge(charge,tol):
    # round up to higher int for .tol
    b,a = math.modf(charge)
    if abs(b)>=tol:
        if b<=0:
            out = a-1
        else:
            out = a+1
    else:
        out = a
    return int(out)

def seperate_fragment(coord,wbo,atcharge):
    nat = coord.shape[0]
    charge0 = round(sum(atcharge))
    if nat ==1:
        # single atom
        fragment_l = [{"coord":coord, "charge":charge0,"atomic_charge":atcharge}]
    else:
        fragment_l = []
        # set check list
        check_l = [False for i in range(nat)]
        for atindex in range(nat):
            if not check_l[atindex]:
                check_l, atomlist = add_fragmentat(atindex+1,wbo,check_l,atomlist=[])
                atomlist = list(set(sorted(atomlist)))
                fcoord = coord.iloc[[i-1 for i in atomlist]]
                fcharge = atcharge.iloc[[i-1 for i in atomlist]]
                fragment_l.append({"coord": fcoord.reset_index(drop=True),
                                   "atomic_charge": fcharge.reset_index(drop=True),
                                   })
        # now set charge based on charge conservation
        charge_floatl = [fragment_l[frindex]["atomic_charge"].loc[:,0].sum() for frindex in range(len(fragment_l))]
        charge_intl = [round_charge(i,0.5) for i in charge_floatl]
        if np.round(sum(charge_intl)) !=np.round(sum(charge_floatl)):
            print(f"please debug here: charge_floatl = {charge_floatl}; charge_intl = {charge_intl}")
        for frindex in range(len(fragment_l)):
            fragment_l[frindex]["charge"] = charge_intl[frindex]
            vnel = get_fr_valenz_nel(fragment_l[frindex]["coord"])
            fragment_l[frindex]["nel"]=vnel+charge_intl[frindex]
    return fragment_l


def creat_randomfolder(wd):
    wd = os.path.abspath(wd)
    random_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
    if not os.path.isdir(os.path.join(wd,random_name)):
        os.makedirs(os.path.join(wd,random_name))
    else:
        rm_folderORfile(os.path.join(wd,random_name))
        os.makedirs(os.path.join(wd,random_name))
    os.chdir(os.path.join(wd,random_name))
    return os.path.join(wd,random_name)



def write_xyz(xyz, name):
    pwd = os.getcwd()
    # xyz.round(10)
    nat = xyz.shape[0]
    with open(f'{pwd}/{name}.xyz', 'w') as f:
        f.write(f'{nat}\n')
        f.write('\n')
        for i in range(nat):
            at, x, y, z = xyz.iloc[i, 0], xyz.iloc[i, 1], xyz.iloc[i, 2], xyz.iloc[i, 3]
            f.write('%2s %12.8f %12.8f %12.8f\n' % (at, x, y, z))
        f.write('\n')


def get_nelectron(xtboutwd):
    with open(f'{xtboutwd}', 'r') as f:
        l = []
        for line in f:
            if "# electrons" in line:
                l.append(line)
        if l:
            e = l[0]
            ne = int(re.findall(r'\-?\d+\.?\d*', e)[0])
        else:
            ne = 0
    return ne


def get_energy(ENERGYFILE):
    speh = np.inf
    if os.path.isfile(f'{ENERGYFILE}'):
        with open(f'{ENERGYFILE}', 'r') as f:
            speh = float(f.readlines()[0])
    return speh



def boltzmanp(a):
    # a: gap or wb(wboe)
    ev2j = 1.602176634 * pow(10, -19)
    return np.exp(-a * ev2j / (1.38 * pow(10, -23) * 1000))


def lowestgeom(infile):
    # input is .log file
    with open(infile, 'r') as f:
        f.seek(0, 0)
        first = f.readline()
        atnr = int(first.split('\n')[0])
        f.seek(0, 0)
        d = {}
        test = f.readlines()
        for i in range(0, len(test), atnr + 2):
            tmp = [test[j] for j in range(i + 1, i + atnr + 2)]
            try:
                energy = float(re.findall(r'\-?\d+\.?\d*', tmp[0].split(':')[1])[0])
            except:
                energy = 'inf'
            d[energy] = tmp
        a = min(d, key=d.get)
        outlist = [str(atnr) + '\n'] + d[a]
        return outlist



def read_topomol(molfile, nat1=None):
    from aRST.script.geom import r0
    # V2000
    with open(molfile) as f:
        context = f.readlines()

    if len(context[3].split())<7:
        nat, nbond = int(''.join(context[3].split()[0][:3])), int(''.join(context[3].split()[0][3:]))
    else:

        nat, nbond = int(context[3].split()[0]), int(context[3].split()[1])
    # coordinates = []
    constrain_l = []
    ele = []

    for i in range(nat):
        # tmp = context[i + 4].split()[:3]
        # coordinates.append([float(k) for k in tmp])
        ele.append(context[i + 4].split()[3])

    for i in range(nbond):
        bond_pair = context[i + nat + 4].split()[:2]
        atom1, atom2 = int(bond_pair[0]), int(bond_pair[1])
        d = r0[ele[atom1 - 1]] + r0[ele[atom2 - 1]]
        tmp_d = {"atoms": {}, "value": {}}
        if nat1 is None:
            tmp_d["atoms"] = [atom1, atom2]
            tmp_d["value"] = d
        else:
            tmp_d["atoms"] = [atom1 + nat1, atom2 + nat1]
            tmp_d["value"] = d
        constrain_l.append(tmp_d)

    return constrain_l, ele

from molbar.barcode import get_molbar_from_file,get_molbar_from_coordinates
def get_molbarcode(coord,molbarlevel=None):
    if coord.shape[0]>1:
        # transform coord to coord_l and ele_l
        coord_l,ele_l = [],[]
        for i in range(coord.shape[0]):
            coord_l.append(coord.loc[i,1:3].tolist())
            ele_l.append(coord.loc[i,0])
        try:
            mb = get_molbar_from_coordinates(coord_l,ele_l)
            if molbarlevel =="loose":
                mb = mb[:4] # only take mb code until Topology Spectrum
        except Exception:
            mb = "unknow"

    else:
        mb = coord.iloc[0,0]
        if mb =="Cl":
            mb = "CL"
    return mb

