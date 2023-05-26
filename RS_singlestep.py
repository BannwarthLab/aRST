import copy
import os
import time
import traceback

import numpy as np
import pandas as pd

"""This script is design to simulate the single-step reactions.
To execute this script:
1. Load module xtb/ml-hessian if the reaction doesn't involve radicals, otherwise you need to pull the xtb version Ying_ml_output from gitlab/Bannwarth
2. Put this script with your already OPTIMIZED reactants geometry files in the same folder.
3. Check if line (15-20) to make sure the calculating parameters are correct (name, charge and multiplicity).
4. Execute the script with command: python3 RS_singlestep.py
"""

################################################
# reactantes setting
reactants = {}
reactants['aa'] = [0, 0]  # reactants name = [charge, uhf]
reactants['b'] = [0, 0]
reaction_type = 0  # here fill in the reaction type, 0: intermolecular reaction; 1: dissociation

tblite_trigger = True
################################################

# covalent radii
r0 = {'H': 0.32, 'D': 0.32, 'He': 0.46, 'Li': 1.20, 'Be': 0.94, 'B': 0.77, 'C': 0.75, 'N': 0.71, 'O': 0.63, 'F': 0.64,
      'Ne': 0.67, 'Na': 1.40, 'Mg': 1.25, 'Al': 1.13, 'Si': 1.04, 'P': 1.10, 'S': 1.02, 'CL': 0.99, 'Ar': 0.96,
      'K': 1.76, 'Ca': 1.54, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11,
      'Ni': 1.10, 'Cu': 1.12, 'Zn': 1.18, 'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17,
      'Rb': 2.10, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
      'Pd': 1.20, 'Ag': 1.28, 'Cd': 1.36, 'In': 1.42, 'Sn': 1.40, 'Sb': 1.40, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31,
      'Cs': 2.32, 'Ba': 1.96, 'La': 1.80, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74, 'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68,
      'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.70, 'Lu': 1.62, 'Hf': 1.52,
      'Ta': 1.46, 'W': 0.95 * 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33, 'Tl': 1.44,
      'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01, 'Ac': 1.86, 'Th': 1.75,
      'Pa': 1.69, 'U': 1.70, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65,
      'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57, 'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34,
      'Mt': 1.29, 'Ds': 1.21, 'Rg': 1.21, 'Cn': 1.22, 'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65,
      'Og': 1.57, 'Cl': 0.99}


# find the active sites
def find_value(dict, fvalue):
    result = False
    for key in dict.keys():
        if dict[key] == fvalue:
            result = True
            break
    return result


def get_key(dict, value):
    return [k for k, v in dict.items() if v == value][0]


def build_gap_d(react_d, h_df, l_df):
    for i in list(h_df.columns):
        for j in list(l_df.columns):
            if i != j:
                struc1 = i + 1
                struc2 = j + 1
                for h_index in h_df.iterrows():
                    for l_index in l_df.iterrows():
                        at1 = h_index[0]
                        at2 = l_index[0]
                        e1 = h_df.loc[at1][i]
                        e2 = l_df.loc[at2][j]
                        gap = e2-e1
                        at1 = at1 +1 # index to atom number
                        at2 = at2 + 1

                        value1 = [f'{struc1},{struc2}', f'{at1},{at2}']
                        value2 = [f'{struc2},{struc1}', f'{at2},{at1}']
                        if gap not in list(react_d.keys()) and not find_value(react_d, value1) and not find_value(
                                react_d, value2) and not np.isnan(gap):
                            react_d[gap] = [f'{struc1},{struc2}', f'{at1},{at2}']
                        elif find_value(react_d, value1) and get_key(react_d, value1) > gap:
                            del react_d[get_key(react_d, value1)]
                            react_d[gap] = [f'{struc1},{struc2}', f'{at1},{at2}']
                        elif find_value(react_d, value2) and get_key(react_d, value2) > gap:
                            del react_d[get_key(react_d, value2)]
                            react_d[gap] = [f'{struc1},{struc2}', f'{at1},{at2}']
    return react_d


def spr_inter_new(reactants, wd_pre, radical_trigger):
    hoaoa_d = {}
    at_l1 = []
    hoaob_d = {}
    at_l2 = []

    luaoa_d = {}
    at_l3 = []
    luaob_d = {}
    at_l4 = []
    for i in range(len(reactants)):
        hoaoa_d[i] = {}
        hoaob_d[i] = {}
        smi = reactants[i]
        if '_' in smi:
            smi = smi.split('_')[0]
        ml_feature = pd.read_csv(f'{wd_pre}/{smi}/ml/ml_feature.csv')
        if 'HOAO_a (eV)' in ml_feature.columns:
            hoaoa_l = list(ml_feature['HOAO_a (eV)'])
            hoaob_l = list(ml_feature['HOAO_b (eV)'])
            target_hoao_a = list(reversed(sorted(hoaoa_l)))[:5]
            target_hoao_b = list(reversed(sorted(hoaob_l)))[:5]
        else:
            hoaoa_l = list(ml_feature['HOAO (eV)'])
            hoaob_l = list(ml_feature['HOAO (eV)'])
            target_hoao_a = list(reversed(sorted(hoaoa_l)))[:5]
            target_hoao_b = list(reversed(sorted(hoaob_l)))[:5]
        for x, y in zip(target_hoao_a, target_hoao_b):
            at_x = hoaoa_l.index(x)
            at_y = hoaob_l.index(y)
            hoaoa_d[i][at_x] = x
            hoaob_d[i][at_y] = y
            at_l1.append(at_x)
            at_l2.append(at_y)
            at_l1 = list(set(at_l1))
            at_l2 = list(set(at_l2))

        luaoa_d[i] = {}
        luaob_d[i] = {}
        if 'LUAO_a (eV)' in ml_feature.columns:
            luaoa_l = list(ml_feature['LUAO_a (eV)'])
            luaob_l = list(ml_feature['LUAO_b (eV)'])
            target_luao_a = list(sorted(luaoa_l))[:5]
            target_luao_b = list(sorted(luaob_l))[:5]
        else:
            luaoa_l = list(ml_feature['LUAO (eV)'])
            luaob_l = list(ml_feature['LUAO (eV)'])
            target_luao_a = list(sorted(luaoa_l))[:5]
            target_luao_b = list(sorted(luaob_l))[:5]
        for x, y in zip(target_luao_a, target_luao_b):
            at_x = luaoa_l.index(x)
            at_y = luaob_l.index(y)
            luaoa_d[i][at_x] = x
            luaob_d[i][at_y] = y
            at_l3.append(at_x)
            at_l4.append(at_y)
            at_l3 = list(set(at_l3))
            at_l4 = list(set(at_l4))
    hoaoa_df = pd.DataFrame(hoaoa_d, index=[i for i in at_l1], columns=[i for i in range(len(reactants))])
    hoaob_df = pd.DataFrame(hoaob_d, index=[i for i in at_l2], columns=[i for i in range(len(reactants))])
    luaoa_df = pd.DataFrame(luaoa_d, index=[i for i in at_l3], columns=[i for i in range(len(reactants))])
    luaob_df = pd.DataFrame(luaob_d, index=[i for i in at_l4], columns=[i for i in range(len(reactants))])

    react_ss = {}
    react_ss = build_gap_d(react_ss, hoaoa_df, luaoa_df)
    react_ss = build_gap_d(react_ss, hoaob_df, luaob_df)

    react_os = {}
    react_os = build_gap_d(react_os, hoaoa_df, luaob_df)
    react_os = build_gap_d(react_os, hoaob_df, luaoa_df)

    if radical_trigger:
        gap_l1 = sorted(list(react_os.keys()))[2:]
        gap_l2 = sorted(list(react_ss.keys()))[2:]
        for u in gap_l1:
            del react_os[u]
        for u in gap_l2:
            del react_ss[u]
        return react_os, react_ss

    else:
        gap_l1 = sorted(list(react_os.keys()))[4:]
        for u in gap_l1:
            del react_os[u]
    return react_os


def spr_intra(smi, wd, disso,tblite):
    # if disso ==ture: dissociate possibility, otherweis intramolecular connection
    # hoao [[h_a],[h_b]]
    # luau [[l_a],[l_b]]
    if '_' in smi:
        smi = smi.split('_')[0]
    if not tblite:
        hl_df = pd.read_csv(f'{wd}/{smi}_ml/ml_feature.csv', index_col=0)
    else:
        hl_df = pd.read_csv(f'{wd}/{smi}_ml/ml_feature_tblite.csv', index_col=0)
    if 'HOAO_a (eV)' in hl_df.columns:
        hoao_l = [hl_df.loc[:, 'HOAO_a (eV)'].tolist(), hl_df.loc[:, 'HOAO_b (eV)'].tolist()]
        luao_l = [hl_df.loc[:, 'LUAO_a (eV)'].tolist(), hl_df.loc[:, 'LUAO_b (eV)'].tolist()]
    else:
        hoao_l = [hl_df.loc[:, 'HOAO (eV)'].tolist(), hl_df.loc[:, 'HOAO (eV)'].tolist()]
        luao_l = [hl_df.loc[:, 'LUAO (eV)'].tolist(), hl_df.loc[:, 'LUAO (eV)'].tolist()]

    # dissociate possibility
    # read wbo of structure and fine the pair with the lowest pair
    # can also returen dis_dirct for the complete information of dissociation
    wbo_df = pd.read_csv(f'{wd}/{smi}_ml/wbo', header=None, sep='\\s+')
    target_pair = [[wbo_df.loc[i, 0], wbo_df.loc[i, 1]] for i in range(wbo_df.shape[0]) if 0.3 < wbo_df.loc[i, 2] < 1.1]
    dis_dirct = {}
    for i in target_pair:
        atom_1 = i[0]
        atom_2 = i[1]
        gap_temp_1 = np.round(abs(hoao_l[0][atom_1 - 1] - luao_l[0][atom_2 - 1]), 2)
        dis_dirct[gap_temp_1] = [atom_1, atom_2]
        dis_dirct[np.round(abs(hoao_l[0][atom_1 - 1] - luao_l[1][atom_2 - 1]), 2)] = [atom_1, atom_2]
        dis_dirct[np.round(abs(hoao_l[1][atom_1 - 1] - luao_l[0][atom_2 - 1]), 2)] = [atom_1, atom_2]
        dis_dirct[np.round(abs(hoao_l[1][atom_1 - 1] - luao_l[1][atom_2 - 1]), 2)] = [atom_1, atom_2]
    dis_gap = min(sorted(dis_dirct.keys()))
    dis_pair = dis_dirct[min(sorted(dis_dirct.keys()))]

    # # connect possibility
    h = copy.deepcopy(hoao_l[0]) + copy.deepcopy(hoao_l[1])
    l = copy.deepcopy(luao_l[0]) + copy.deepcopy(luao_l[1])
    h.sort()
    l.sort()
    i, j, gap = 0, 0, float('inf')
    while i < len(h) and j < len(l):
        gap = min(gap, abs(h[i] - l[j]))
        a_1 = i
        a_2 = j
        if h[i] > l[j]:
            j += 1
        else:
            i += 1
        if h[a_1] in hoao_l[0]:
            atom_a = hoao_l[0].index(h[a_1])
        else:
            atom_a = hoao_l[1].index(h[a_1])
        if l[a_2] in luao_l[0]:
            atom_b = luao_l[0].index(l[a_2])
        else:
            atom_b = luao_l[1].index(l[a_2])
        if not disso:
            return gap, [atom_a, atom_b]
        else:
            return dis_gap, dis_pair


class Coord:
    def __init__(self, coord):
        self.coord = coord

    def get_principal_moment_vec(self, at):
        t1, t2, t3, t4, t5, t6 = 0, 0, 0, 0, 0, 0
        # matrix for moments of inertia is of form, with CN in mass as weight
        # *           |   y**2+z**2                         |
        # *           |    -y*x       z**2+x**2             | -i =0
        # *           |    -z*x        -z*y       x**2+y**2 |
        xv, yv, zv = 0, 0, 0
        for i in range(self.coord.shape[0]):
            if i != at - 1:
                x = self.coord.iloc[i, 1] - self.coord.iloc[at - 1, 1]
                y = self.coord.iloc[i, 2] - self.coord.iloc[at - 1, 2]
                z = self.coord.iloc[i, 3] - self.coord.iloc[at - 1, 3]
                r = np.linalg.norm(np.array([x, y, z]))
                rcov_c = float(r0[f'{self.coord.iloc[at - 1, 0]}'])
                rcov_i = float(r0[f'{self.coord.iloc[i, 0]}'])
                rr = (rcov_i + rcov_c) / r
                cn = 1 / (1 + np.exp(- 16 * rr - 1))
                atmass = 0.5 * cn / (r ** 2)
                t1 = t1 + atmass * (y ** 2 + z ** 2)
                t2 = t2 + atmass * (-y * x)
                t3 = t3 + atmass * (z ** 2 + x ** 2)
                t4 = t4 + atmass * (-z * x)
                t5 = t5 + atmass * (-z * y)
                t6 = t6 + atmass * (x ** 2 + y ** 2)
                xv = xv + atmass * x
                yv = yv + atmass * y
                zv = zv + atmass * z
        # print(cn)
        mv = np.array([xv, yv, zv])
        arep = np.zeros((3, 3))
        arep[0][0], arep[1][0], arep[1][1], arep[2][0], arep[2][1], arep[2][2] = t1, t2, t3, t4, t5, t6
        arep = arep + arep.T - np.diag(arep.diagonal())

        v, P = np.linalg.eigh(arep)
        #     print(v.round(3))
        #     print(P.round(3))
        return mv, P

    def get_punkt_coord(self, at):
        a = np.array([self.coord.iloc[at - 1, 1],
                      self.coord.iloc[at - 1, 2],
                      self.coord.iloc[at - 1, 3]])
        return a

    def get_nat(self):
        return self.coord.shape[0]

    def test_if_right_handed_system(self, at):
        P = Coord.get_principal_moment_vec(self, at)[1]
        a, b, c = P[:, 0], P[:, 1], P[:, 2]
        cross_p = [np.round(j, 3) for j in np.cross(b, c)]
        dot_p = np.dot(a, cross_p)
        if dot_p > 0:
            return True
        else:
            return False

    def geom_orientation(self, at):
        P = Coord.get_principal_moment_vec(self, at)[1]
        mv = Coord.get_principal_moment_vec(self, at)[0]
        v1 = np.array(P[:, 2])  # pz
        v2 = mv
        r1 = np.linalg.norm(np.array(v1))
        r2 = np.linalg.norm(np.array(v2))
        dot_p = np.dot(v1, v2)
        c = dot_p / (r1 * r2)
        if c >= 0:
            return True
        else:
            return False

    def shift(self, at, d):
        # shift the molecule with aktiv as center to (0,0,d)
        x = self.coord.iloc[at - 1, 1]
        y = self.coord.iloc[at - 1, 2]
        z = self.coord.iloc[at - 1, 3] - d
        for i in range(self.coord.shape[0]):
            self.coord.iloc[i, 1] = self.coord.iloc[i, 1] - x
            self.coord.iloc[i, 2] = self.coord.iloc[i, 2] - y
            self.coord.iloc[i, 3] = self.coord.iloc[i, 3] - z
        return self.coord

    def rot_y(self):
        # right-handed cs and rot
        beta = np.deg2rad(180)
        Ry = np.array([[np.cos(beta), 0, np.sin(beta)],
                       [0, 1, 0],
                       [-np.sin(beta), 0, np.cos(beta)]])
        for i in range(self.coord.shape[0]):
            x, y, z = self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3]
            a = [x, y, z]
            a_new = np.dot(Ry, a)
            self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3] = a_new[0], a_new[1], a_new[2]
        return self.coord

    def rot_x(self):
        # right-handed cs and rot
        alpha = np.deg2rad(180)
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(alpha), -np.sin(alpha)],
                       [0, np.sin(alpha), np.cos(alpha)]])
        for i in range(self.coord.shape[0]):
            x, y, z = self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3]
            a = [x, y, z]
            a_new = np.dot(Rx, a)
            self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3] = a_new[0], a_new[1], a_new[2]
        return self.coord

    def rot_z(self):
        # right-handed cs and rot
        gamma = np.deg2rad(180)
        Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                       [np.sin(gamma), np.cos(gamma), 0],
                       [0, 0, 1]])
        for i in range(self.coord.shape[0]):
            x, y, z = self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3]
            a = [x, y, z]
            a_new = np.dot(Rz, a)
            self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3] = a_new[0], a_new[1], a_new[2]
        return self.coord

    def change_cs(self, P):
        for i in range(self.coord.shape[0]):
            x, y, z = self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3]
            a = [x, y, z]
            a_new = np.dot(P, a)
            self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3] = a_new[0], a_new[1], a_new[2]
        return self.coord

    def l2r(self):
        for i in range(self.coord.shape[0]):
            x, y, z = self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3]
            a = [x, y, z]
            a_new = [x, -y, z]
            self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3] = a_new[0], a_new[1], a_new[2]
        return self.coord


def alignment(coord1, coord2, at1, at2):
    C1 = Coord(coord1)
    P1 = C1.get_principal_moment_vec(at=at1)[1]
    rhanded = C1.test_if_right_handed_system(at=at1)
    orient = C1.geom_orientation(at=at1)
    C1 = C1.shift(at=at1, d=0)
    if not rhanded:
        C1 = Coord(C1).l2r()  # flip y to -y
    if orient:
        C1 = Coord(C1).rot_x()  # rot around x axis 180
    C1 = Coord(C1).change_cs(P1)
    name1 = C1.iloc[at1 - 1, 0]

    C2 = Coord(coord2)
    P2 = C2.get_principal_moment_vec(at=at2)[1]
    rhanded = C2.test_if_right_handed_system(at=at2)
    orient = C2.geom_orientation(at=at2)
    C2 = C2.shift(at=at2, d=0)
    if not rhanded:
        C2 = Coord(C2).l2r()  # flip y to -y
    if orient:
        C2 = Coord(C2).rot_x()  # rot around x axis 180
    C2 = Coord(C2).change_cs(P2)
    name2 = C2.iloc[at2 - 1, 0]
    rr = (r0[name1] + r0[name2])
    d = np.round(rr * 5, 2)
    C2 = Coord(C2).shift(at=at2, d=d)

    coord_merge = pd.concat([C1, C2])
    return coord_merge, rr


def creat_folder(wd):
    if not os.path.isdir(f'{wd}'):
        os.makedirs(f'{wd}')
    os.chdir(f'{wd}')
    return


def write_xyz(xyz, name):
    xyz.to_csv(f'{name}.xyz', sep='\t', index=False, header=None)
    nat = xyz.shape[0]
    with open(f'{name}.xyz', 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(str(nat) + '\n' + '\n' + content)





def sel_xtbcoord(wd, react):
    if '_' in react:
        react = react.split('_')[0]
    if os.path.isfile(f'{wd}/{react}_ml/xtbopt.xyz'):
        xyz = pd.read_table(f'{wd}/{react}_ml/xtbopt.xyz', header=None, skiprows=2, sep='\s+')
    else:
        if os.path.isfile(f'{wd}/{react}_ml/xtblast.xyz'):
            xyz = pd.read_table(f'{wd}/{react}_ml/xtbopt.xyz', header=None, skiprows=2, sep='\s+')
        else:
            xyz = pd.read_table(f'{wd}/{react}_ml/coord.xyz', header=None, skiprows=2, sep='\s+')
    return xyz


def write_scan_input(at1, at2, d1, d2):
    pwd = os.getcwd()
    with open('scan.inp', 'w') as f3:
        f3.write(f'$constrain\n'
                 f'force constant=1\n'
                 f'$scan\n'
                 f'distance: {at1},{at2},{d1};{d1},{d2},30\n$end')


def push(push_dict, react_l, wd, wd_2):
    # input push_dict,reactants_list, working directory(i use it to search xyz file),wd_2: directory you want to place the merged xyz and scan.inp
    struc_1, struc_2 = push_dict[i][0].split(',')
    struc_1, struc_2 = react_l[int(struc_1) - 1], react_l[int(struc_2) - 1]
    at_1, at_2 = push_dict[i][1].split(',')
    at_1, at_2 = int(at_1), int(at_2)

    coord_1 = sel_xtbcoord(wd, struc_1)
    coord_2 = sel_xtbcoord(wd, struc_2)
    nat_1 = Coord(coord_1).get_nat()

    coord_merge, rr = alignment(coord_1, coord_2, at_1, at_2)
    nreact = str(count)
    path_l.append(nreact)
    creat_folder(f'{wd_2}')
    write_xyz(coord_merge, 'coord')
    charg = int(rec[struc_1]['charge'] + rec[struc_2]['charge'])
    uhf = int(rec[struc_1]['uhf'] + rec[struc_2]['uhf'])
    write_scan_input(at_1, at_2 + nat_1, rr * 5, rr)


def pull(dis_dict, wd, wd_2):
    creat_folder(f'{wd_2}')
    mdgap = min(dis_dict.keys())
    react_name = dis_dict[mdgap][0]
    dis_pair = dis_dict[mdgap][1]
    coord = sel_xtbcoord(wd, react_name)
    r1 = Coord(coord).get_punkt_coord(dis_pair[0])
    r2 = Coord(coord).get_punkt_coord(dis_pair[1])
    r = np.linalg.norm(r1 - r2)
    write_xyz(coord, 'coord')
    write_scan_input(dis_pair[0], dis_pair[1], r, r * 5)


def ml_calculation(tblite_trigger,charg, uhf):
    if tblite_trigger:
        os.system(f'tblite run coord.xyz --xtbml_xyz -c {charg} --spin {uhf} > xtb.out')
    else:
        os.system(f'xtb coord.xyz --ml_feature -c {charg} --uhf {uhf} > xtb.out')
    time.sleep(1)
    return





########################################
rec = {}
wd = '/home/guests/ychen/test/test_sp'
# recording reactants info
for i in range(len(reactants.keys())):
    smi = list(reactants.keys())[i][0]
    if type(smi) ==tuple: # duplicated reactants
        for j in range(int(smi[1])):
            smi = smi[0]
            smi_tmp = smi + f'_{j}'
            rec[smi_tmp] = {}
            rec[smi_tmp]['charge'] = reactants[list(reactants.keys())[i]][0]
            rec[smi_tmp]['uhf'] = reactants[list(reactants.keys())[i]][1]
    else:
        smi = list(reactants.keys())[i]
        rec[smi] = {}
        rec[smi]['charge'] = reactants[smi][0]
        rec[smi]['uhf'] = reactants[smi][1]

react_l = list(rec.keys())

# ml features
for name in react_l:
    if '_' in name:
        name_tmp = name.split('_')[0]
        C = Coord(pd.read_table(f'{wd}/{name_tmp}.xyz', header=None, skiprows=2, sep='\s+'))
    else:
        C = Coord(pd.read_table(f'{wd}/{name}.xyz', header=None, skiprows=2, sep='\s+'))
    creat_folder(f'{wd}/{name}_ml')
    write_xyz(C.coord, 'coord')
    charg = rec[name]['charge']
    uhf = rec[name]['uhf']
    ml_calculation(tblite_trigger,charg,uhf)
    while not os.path.isfile(f'{wd}/{name}_ml/xtb.out'):
        time.sleep(5)
    else:
        C.coord.columns = ['at', 'x', 'y', 'z']
        if tblite_trigger:
            ml_feature = pd.read_csv(f'{wd}/{name}_ml/ml_feature_tblite.csv')
        else:
            ml_feature = pd.read_csv(f'{wd}/{name}_ml/ml_feature.csv')
        rec[name]['coord'] = C.coord

        # at_features = pd.concat([C.coord, c, cn], axis=1)


if reaction_type == 0:
    ## inter possibility (no radicals: top 4; radicals: top 2)
    try:
        radical_trigger = False
        for name in react_l:
            if rec[name]['uhf'] != 0:
                radical_trigger = True
        if not radical_trigger:
            push_dict_1 = spr_inter_new(reactants=react_l, wd_pre=f'{wd}', radical_trigger=False,tblite=tblite_trigger)
            # {11.11141: ['2,1', '1,2'], 11.3398: ['2,1', '1,1'], 9.774439999999998: ['2,1', '2,2'], 10.00283: ['2,1', '2,1'], 11.00347: ['1,2', '2,1']}
            gap_l = sorted(list(push_dict_1.keys()))[:4]
        else:
            push_dict_1, push_dict_2 = spr_inter_new(reactants=react_l, wd_pre=f'{wd}', radical_trigger=True,tblite=tblite_trigger)
            gap_l = sorted(list(push_dict_1.keys()))[:2] + sorted(list(push_dict_2.keys()))[:2]
            push_dict_1.update(push_dict_2)

        path_l = []
        count = 1
        for i in gap_l:
            # get reactants name and active atoms
            struc_1, struc_2 = push_dict_1[i][0].split(',')
            struc_1, struc_2 = react_l[int(struc_1) - 1], react_l[int(struc_2) - 1]
            at_1, at_2 = push_dict_1[i][1].split(',')
            at_1, at_2 = int(at_1), int(at_2)

            coord_1 = (sel_xtbcoord(wd, struc_1))
            coord_2 = (sel_xtbcoord(wd, struc_2))
            nat_1 = Coord(coord_1).get_nat()
            coord_merge, rr = alignment(coord_1, coord_2, at_1, at_2)
            nreact = str(count)
            path_l.append(nreact)
            creat_folder(f'{wd}/{nreact}')
            write_xyz(coord_merge, 'coord')
            charg = int(rec[struc_1]['charge'] + rec[struc_2]['charge'])
            uhf = int(rec[struc_1]['uhf'] + rec[struc_2]['uhf'])
            write_scan_input(at_1, at_2 + nat_1, rr * 5, rr)
            count += 1
    except:
        print(traceback.print_exc())
        pass
else:
    if reaction_type == 1:
        dis_dict = {}
        for i in react_l:
            dis_gap, dis_pair = spr_intra(smi=i, wd=wd, disso=True,tblite=tblite_trigger)
            dis_dict[dis_gap] = i, dis_pair

        creat_folder(f'{wd}/scan')
        mdgap = min(dis_dict.keys())
        react_name = dis_dict[mdgap][0]
        dis_pair = dis_dict[mdgap][1]
        coord = sel_xtbcoord(wd, react_name)
        r1 = Coord(coord).get_punkt_coord(dis_pair[0])
        r2 = Coord(coord).get_punkt_coord(dis_pair[1])
        r = np.linalg.norm(r1 - r2)
        write_xyz(coord, 'coord')
        write_scan_input(dis_pair[0], dis_pair[1], r, r * 5)
    else:
        print("reaction type wrong")

