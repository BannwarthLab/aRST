import copy
import os
import shutil

import numpy as np
import pandas as pd

"""This script is design to simulate the single-step reactions.
To execute this script:
1. Load module xtb/ml-hessian if the reaction doesn't involve radicals, otherwise you need to pull the xtb version Ying_ml_output from gitlab/Bannwarth
2. Put this script with your already OPTIMIZED reactants geometry files in the same folder.
3. Check if line (15-20) to make sure the calculating parameters are correct (name, charge and multiplicity).
4. Execute the script with command: python3 RS_singlestep.py
"""

####################### Reaction parameters #########################

#reaction_type = 1  # here fill in the reaction type, 0: intermolecular reaction; 1: dissociation

####################### Functions #########################

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
def spr_inter(reactants, wd):
    # 1: hexyne	2: H2
    # {'2,1': [(9.61370222, [1, 15]), (9.73839884, [1, 14]), (15.28320733, [1, 13])]}
    if len(reactants) > 1:
        h = []
        l = []
        for smi in reactants:
            df = pd.read_csv(f'{wd}/{smi}_ml/ml_feature.csv')
            if 'HOAO_a (eV)' in df.columns:
                h.append(df['HOAO_a (eV)'].tolist())
                h.append(df['HOAO_b (eV)'].tolist())
                l.append(df['LUAO_a (eV)'].tolist())
                l.append(df['LUAO_b (eV)'].tolist())
            else:
                h.append(df['HOAO (eV)'].tolist())
                h.append(df['HOAO (eV)'].tolist())
                l.append(df['LUAO (eV)'].tolist())
                l.append(df['LUAO (eV)'].tolist())
        d_ss = {}
        d_os = {}
        for i in range(len(h)):
            for j in range(len(l)):
                struc_1_nr = int(i / 2) + 1
                struc_2_nr = int(j / 2) + 1
                if struc_1_nr != struc_2_nr:
                    h_i = sorted(h[i])[:5]
                    l_j = sorted(l[j])[:5]

                    # Setting up the iterative dictionary
                    if f'{struc_1_nr},{struc_2_nr}' not in d_ss.keys():
                        d_ss[f'{struc_1_nr},{struc_2_nr}'] = {}
                    if f'{struc_1_nr},{struc_2_nr}' not in d_os.keys():
                        d_os[f'{struc_1_nr},{struc_2_nr}'] = {}

                    for one in range(len(h_i)):
                        for two in range(len(l_j)):
                            # Exclusion of hydrogen ions
                            try:
                                h_i[one] = float(h_i[one])
                                l_j[two] = float(l_j[two])
                            except ValueError:
                                h_i[one] = float('inf')
                                l_j[two] = float('inf')

                            key = np.round(float(l_j[two])-float(h_i[one]), 3)
                            at_nr_origial_1 = h[i].index(h_i[one]) + 1
                            at_nr_origial_2 = l[j].index(l_j[two]) + 1

                            if (i % 2 == 0 and j % 2 == 0) or (i % 2 != 0 and j % 2 != 0):
                                if key in d_ss[f'{struc_1_nr},{struc_2_nr}'].keys() \
                                        and [at_nr_origial_1, at_nr_origial_2] not in d_ss[f'{struc_1_nr},{struc_2_nr}'][key]:
                                    d_ss[f'{struc_1_nr},{struc_2_nr}'][key].append([at_nr_origial_1, at_nr_origial_2])
                                else:
                                    d_ss[f'{struc_1_nr},{struc_2_nr}'][key] = [[at_nr_origial_1, at_nr_origial_2]]
                            else:
                                if key in d_os[f'{struc_1_nr},{struc_2_nr}'].keys() \
                                        and [at_nr_origial_1, at_nr_origial_2] not in d_os[f'{struc_1_nr},{struc_2_nr}'][key]:
                                    d_os[f'{struc_1_nr},{struc_2_nr}'][key].append([at_nr_origial_1, at_nr_origial_2])
                                else:
                                    d_os[f'{struc_1_nr},{struc_2_nr}'][key] = [[at_nr_origial_1, at_nr_origial_2]]
                            # First filtering of degenerate pairs
                            for g_1 in list(d_ss[f'{struc_1_nr},{struc_2_nr}'].keys()):
                                if 0 < abs(key - g_1) < 0.01:
                                    temp_pair = list(d_ss[f'{struc_1_nr},{struc_2_nr}'][g_1])
                                    temp_pair[0].append('*')
                                    d_ss[f'{struc_1_nr},{struc_2_nr}'][key] = d_ss[f'{struc_1_nr},{struc_2_nr}'][key] + temp_pair
                                    del d_ss[f'{struc_1_nr},{struc_2_nr}'][g_1]

                            for g_2 in list(d_os[f'{struc_1_nr},{struc_2_nr}'].keys()):
                                if 0 < abs(key - g_2) < 0.01:
                                    temp_pair = list(d_os[f'{struc_1_nr},{struc_2_nr}'][g_2])
                                    temp_pair[0].append('*')
                                    d_os[f'{struc_1_nr},{struc_2_nr}'][key] = d_os[f'{struc_1_nr},{struc_2_nr}'][key] + temp_pair
                                    del d_os[f'{struc_1_nr},{struc_2_nr}'][g_2]
        # print(d_ss)
        # print(d_os)
        # delete excuss data
        for k_1 in list(d_ss.keys()):
            gap_list = sorted(list(d_ss[k_1].keys()))[2:]
            for l_1 in gap_list:
                del d_ss[k_1][l_1]
        for k_2 in list(d_os.keys()):
            gap_list = sorted(list(d_os[k_2].keys()))[2:]
            for l_2 in gap_list:
                del d_os[k_2][l_2]

    print(d_ss)
    return d_ss,d_os


def spr_intra(smi, wd, disso):
    # if disso ==ture: dissociate possibility, otherweis intramolecular connection
    # hoao [[h_a],[h_b]]
    # luau [[l_a],[l_b]]

    hl_df = pd.read_csv(f'{wd}/{smi}_ml/ml_feature.csv', index_col=0)
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
        print(f'{smi}: {dis_gap}, {dis_pair}')
        return dis_gap, dis_pair


def shift(df, aktiv, d):
    # shift the molecule with aktiv as center to (0,0,d)
    x = df.iloc[aktiv, 1]
    y = df.iloc[aktiv, 2]
    z = df.iloc[aktiv, 3] - d
    for i in range(df.shape[0]):
        df.iloc[i, 1] = df.iloc[i, 1] - x
        df.iloc[i, 2] = df.iloc[i, 2] - y
        df.iloc[i, 3] = df.iloc[i, 3] - z
    return df

def flip(df, aktiv):
    x = df.iloc[aktiv, 1]
    y = df.iloc[aktiv, 2]
    z = df.iloc[aktiv, 3] - d
    for i in range(df.shape[0]):
        df.iloc[i, 1] = - x
        df.iloc[i, 2] = - y
        df.iloc[i, 3] = - z
    return df



def get_principal_moment(coord, aktiv):
    r0 = {'H': 0.32, 'D': 0.32, 'He': 0.46, 'Li': 1.20, 'Be': 0.94, 'B': 0.77, 'C': 0.75, 'N': 0.71, 'O': 0.63,
          'F': 0.64,
          'Ne': 0.67, 'Na': 1.40, 'Mg': 1.25, 'Al': 1.13, 'Si': 1.04, 'P': 1.10, 'S': 1.02, 'CL': 0.99, 'Ar': 0.96,
          'K': 1.76, 'Ca': 1.54, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11,
          'Ni': 1.10, 'Cu': 1.12, 'Zn': 1.18, 'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17,
          'Rb': 2.10, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
          'Pd': 1.20, 'Ag': 1.28, 'Cd': 1.36, 'In': 1.42, 'Sn': 1.40, 'Sb': 1.40, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31,
          'Cs': 2.32, 'Ba': 1.96, 'La': 1.80, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74, 'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68,
          'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.70, 'Lu': 1.62, 'Hf': 1.52,
          'Ta': 1.46, 'W': 0.95 * 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33,
          'Tl': 1.44,
          'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01, 'Ac': 1.86, 'Th': 1.75,
          'Pa': 1.69, 'U': 1.70, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65,
          'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57, 'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34,
          'Mt': 1.29, 'Ds': 1.21, 'Rg': 1.21, 'Cn': 1.22, 'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65,
          'Og': 1.57, 'Cl': 0.99}
    t1, t2, t3, t4, t5, t6 = 0, 0, 0, 0, 0, 0
    # matrix for moments of inertia is of form, with CN in mass as weight
    # *           |   y**2+z**2                         |
    # *           |    -y*x       z**2+x**2             | -i =0
    # *           |    -z*x        -z*y       x**2+y**2 |
    cn = 0
    for i in range(coord.shape[0]):
        if i != aktiv:
            x = coord.iloc[i, 1] - coord.iloc[aktiv, 1]
            y = coord.iloc[i, 2] - coord.iloc[aktiv, 2]
            z = coord.iloc[i, 3] - coord.iloc[aktiv, 3]
            r = np.linalg.norm(np.array([x, y, z]))
            rcov_c = float(r0[f'{coord[0][aktiv]}'])
            rcov_i = float(r0[f'{coord[0][i]}'])
            rr = (rcov_i + rcov_c) / r
            cn = 1 / (1 + np.exp(- 16 * (rr - 1)))
            atmass = 0.5 * cn / (r ** 2)
            t1 = t1 + atmass * (y ** 2 + z ** 2)
            t2 = t2 + atmass * (-y * x)
            t3 = t3 + atmass * (z ** 2 + x ** 2)
            t4 = t4 + atmass * (-z * x)
            t5 = t5 + atmass * (-z * y)
            t6 = t6 + atmass * (x ** 2 + y ** 2)
    # print(cn)
    arep = np.zeros((3, 3))
    arep[0][0], arep[1][0], arep[1][1], arep[2][0], arep[2][1], arep[2][2] = t1, t2, t3, t4, t5, t6
    arep = arep + arep.T - np.diag(arep.diagonal())

    v, P = np.linalg.eigh(arep)
    #     print(v.round(3))
    #     print(P.round(3))
    return v, P


def test_if_right_handed_system(P):
    cross_p1 = [np.round(j, 3) for j in np.cross(P[0], P[1])]
    cross_p2 = [np.round(j, 3) for j in np.cross(P[1], P[0])]
    P_z = [np.round(j, 3) for j in P[2]]
    #     print(cross_p1)
    #     print(cross_p2)
    #     print(P[2])

    i = 0
    while i < len(cross_p1):
        if P_z[i] != cross_p1[i]:
            if cross_p2[i] == P_z[i]:
                trigger = True
                i = i + 1
                pass
            else:
                print("error in main axes")
                break
        else:
            trigger = False
            i = i + 1
    if trigger:
        return True
    else:
        return False


def transform(coord, P):
    for i in range(coord.shape[0]):
        x, y, z = coord.iloc[i, 1], coord.iloc[i, 2], coord.iloc[i, 3]
        a = [x, y, z]
        a_new = np.dot(a, P)
        coord.iloc[i, 1], coord.iloc[i, 2], coord.iloc[i, 3] = a_new[0], a_new[1], a_new[2]
    return coord


def aligment(smi_a, smi_b, rec, d, wd, trick):
    """"Diagonalized inertia tensor: P^-1 I P = I'
            rotate the I: P_new^-1 IQ P_new =I'
            since P_new = unit matrix
            Q = P
            Equivalent to use P as rotation matrix"""
    # alignment for first molecule

    # read coord date
    if os.path.isfile(f'{wd}/{smi_a}_ml/xtbopt.xyz'):
        coord_a = pd.read_table(f'{wd}/{smi_a}_ml/xtbopt.xyz', header=None, skiprows=2, sep='\s+')
    else:
        coord_a = pd.read_table(f'{wd}/{smi_a}_ml/{smi_a}.xyz', header=None, skiprows=2, sep='\s+')
    coord_a = coord_a.round(3)

    # get local moment of inertia and get the eigen vector matrix
    v1, P1 = get_principal_moment(coord=coord_a, aktiv=rec[0] - 1)
    # test if eigen vector matrix is right handed
    if test_if_right_handed_system(P1):
        coord_temp = transform(coord_a.copy(deep=True), P1)
        vt, Pt = get_principal_moment(coord=coord_temp, aktiv=rec[0] - 1)
        P1 = np.dot(P1, Pt.T)
    coord_a = transform(coord_a, P1)
    v1, P1 = get_principal_moment(coord=coord_a, aktiv=rec[0] - 1)
    coord_a = shift(coord_a, rec[0] - 1, 0)

    # alignment for second molecule
    if os.path.isfile(f'{wd}/{smi_b}_ml/xtbopt.xyz'):
        coord_b = pd.read_table(f'{wd}/{smi_b}_ml/xtbopt.xyz', header=None, skiprows=2, sep='\s+')
    else:
        coord_b = pd.read_table(f'{wd}/{smi_b}_ml/{smi_b}.xyz', header=None, skiprows=2, sep='\s+')
    coord_b = coord_b.round(3)
    v2, P2 = get_principal_moment(coord=coord_b, aktiv=rec[1] - 1)
    # test if eigen vector matrix is right handed
    if test_if_right_handed_system(P2):
        coord_temp = transform(coord_b.copy(deep=True), P2)
        vt, Pt = get_principal_moment(coord=coord_temp, aktiv=rec[1] - 1)
        P2 = np.dot(P2, Pt.T)

    # largest main axis point into negativ
    Q = np.zeros((3, 3))
    Q[0][0], Q[1][1], Q[2][2] = 1, 1, -1
    P2 = np.dot(P2, Q.T)

    coord_b = transform(coord_b, P2)
    v2, P2 = get_principal_moment(coord=coord_b, aktiv=rec[1] - 1)
    if trick:
        coord_b = flip(coord_b,rec[1] - 1)
    coord_b = shift(coord_b, rec[1] - 1, d)

    # merge two coord files

    product = pd.concat([coord_a, coord_b])
    atomnr_1 = coord_a.shape[0]
    atmonr_2 = coord_b.shape[0]

    atomnr = str(atomnr_1 + atmonr_2)
    product.to_csv('product.xyz', sep='\t', index=False, header=None)
    with open('product.xyz', 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(atomnr + '\n' + '\n' + content)
    with open('cadm.bild', 'w+') as f2:
        f2.write('.color red\n')
        f2.write(f'.arrow 0 0 0 0 0 10 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 0 {P1[2][0] * 4} {P1[2][1] * 4} {P1[2][2] * 4} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 0 {P1[0][0]} {P1[0][1]} {P1[0][2]} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 0 {P1[1][0] * 2} {P1[1][1] * 2} {P1[1][2] * 2} 0.05 0.150000\n')

        f2.write('.color green\n')
        f2.write(f'.arrow 0 0 0 10 0 0 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[2][0] * 4} {P2[2][1] * 4} {(P2[2][2] * 4 + d)} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[0][0]} {P2[0][1]} {P2[0][2] + d} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[1][0] * 2} {P2[1][1] * 2} {(P2[1][2] * 2 + d)} 0.05 0.150000\n')

        f2.write('.color blue\n')
        f2.write(f'.arrow 0 0 0 0 10 0 0.05 0.150000\n')
    return atomnr_1


def push(reactants_temp, push_rec, push_smi_a, push_smi_b, wd):
    # get reactant name
    push_smi_a = reactants_temp[int(push_smi_a) - 1]
    push_smi_b = reactants_temp[int(push_smi_b) - 1]
    # creat folder for aligned molecule
    if not os.path.isdir(f'{wd}/scan'):
        os.makedirs(f'{wd}/scan')
    os.chdir(f'{wd}/scan')

    # get covalent bond length of the active sites, relativ distance between active sites is 5*r_colvalent
    coord_1 = pd.read_table(f'{wd}/{push_smi_a}.xyz', header=None, skiprows=2, sep='\s+')
    coord_2 = pd.read_table(f'{wd}/{push_smi_b}.xyz', header=None, skiprows=2, sep='\s+')
    r_col_1 = r0[coord_1[0][int(push_rec[0][0]) - 1]]
    r_col_2 = r0[coord_2[0][int(push_rec[0][1]) - 1]]
    align_distance = np.round((r_col_1 + r_col_2) * 5, 2)

    # generate aligned molecules and return atomic number of first molecule
    atomnnr_a = aligment(smi_a=push_smi_a, smi_b=push_smi_b, rec=push_rec[0], d=align_distance, wd=wd,trick=False)
    if len(push_rec) > 1:
        # concerterted mode
        push_rec_p = [[push_rec[0][0], push_rec[0][1] + atomnnr_a], [push_rec[1][0], push_rec[1][1] + atomnnr_a]]
        start_distance_1 = align_distance
        end_distance_1 = np.round((r_col_1 + r_col_2) * 0.9, 2)

        # get relativ distance/covalent bond length of second pairs
        coord = pd.read_table(f'{wd}/scan/product.xyz', header=None, skiprows=2, sep='\s+')
        x = coord[1][push_rec[1][0] - 1] - coord[1][push_rec[1][1] - 1 + atomnnr_a]
        y = coord[2][push_rec[1][0] - 1] - coord[2][push_rec[1][1] - 1 + atomnnr_a]
        z = coord[3][push_rec[1][0] - 1] - coord[3][push_rec[1][1] - 1 + atomnnr_a]
        start_distance_2 = np.round(np.sqrt(x ** 2 + y ** 2 + z ** 2) * 1.1, 2)
        r_col_3 = r0[coord_1[0][int(push_rec[1][0]) - 1]]
        r_col_4 = r0[coord_2[0][int(push_rec[1][1]) - 1]]
        end_distance_2 = np.round((r_col_3 + r_col_4) * 0.9, 2)
        distance_rec = [start_distance_1, end_distance_1, start_distance_2, end_distance_2]
    else:
        # default mode
        push_rec_p = [[push_rec[0][0], push_rec[0][1] + atomnnr_a]]
        start_distance = align_distance
        end_distance = np.round((r_col_1 + r_col_2) * 0.9, 2)
        distance_rec = [start_distance, end_distance]

    push_process(push_rec_p, distance_rec)
    return


def push_process(push_rec_p, distance_rec):
    if len(push_rec_p) > 1:
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant=1\n'
                     f'$scan\n'
                     f'mode = concerted\n'
                     f'distance: {push_rec_p[0][0]},{push_rec_p[0][1]},{distance_rec[0]};{distance_rec[0]},{distance_rec[1]},30\n'
                     f'distance: {push_rec_p[1][0]},{push_rec_p[1][1]},{distance_rec[2]};{distance_rec[2]},{distance_rec[3]},30\n$end')
            # f'distance: {push_rec_p[1][0]},{push_rec_p[1][1]},{distance_rec[2]};{distance_rec[0]},{distance_rec[3]},30\n$end')
    else:
        with open('scan.inp', 'w') as f3:
            f3.write(f'$constrain\n'
                     f'force constant=1\n'
                     f'$scan\n'
                     f'distance: {push_rec_p[0][0]},{push_rec_p[0][1]},{distance_rec[0]};{distance_rec[0]},{distance_rec[1]},30\n$end')
    return


def atom_distance(coord, nr_1, nr_2):
    x = coord[1][int(nr_1) - 1] - coord[1][int(nr_2) - 1]
    y = coord[2][int(nr_1) - 1] - coord[2][int(nr_2) - 1]
    z = coord[3][int(nr_1) - 1] - coord[3][int(nr_2) - 1]
    distance = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    # if max([coord_1[1][int(nr_1)-1], coord_1[2][int(nr_1)-1], coord_1[3][int(nr_1)-1]]) < 1:
    #      distance = distance * 10

    return distance


def pull(reactants_temp, pull_rec, pull_smi_a, wd):
    # get reactant name
    pull_smi_a = reactants_temp[int(pull_smi_a) - 1]
    # creat folder for aligned molecule
    if not os.path.isdir(f'{wd}/scan'):
        os.makedirs(f'{wd}/scan')
    os.chdir(f'{wd}/scan')

    # get current distance of the active sites
    coord = pd.read_table(f'{wd}/{pull_smi_a}.xyz', header=None, skiprows=2, sep='\s+')
    start_distance = atom_distance(coord,pull_rec[0],pull_rec[1])
    end_distance = start_distance * 5

    # generate scan files
    shutil.copy(f'{wd}/{pull_smi_a}.xyz',f'{wd}/scan')
    with open('scan.inp', 'w') as f3:
        f3.write(f'$constrain\n'
                 f'force constant=2\n'
                 f'$scan\n'
                 f'distance: {pull_rec[0]},{pull_rec[1]},{start_distance};{start_distance},{end_distance},30\n$end')
    return


#########################################################


def call_dissociation(target_name_1:str, c_1:int, u_1:int):
    wd = os.getcwd()
    d = {}
    d[target_name_1] = [c_1, u_1]
    #if target_name_2:
        # d[target_name_2] = [c_2, u_2]

    #call xTB
    #os.system('module load xtb/ml-hessian')

    for i in list(d.keys()):
        if os.path.isdir(f'{wd}/{i}_ml'):
                shutil.rmtree(f'{wd}/{i}_ml')
        os.makedirs(f'{wd}/{i}_ml')
        os.system(f'cp {wd}/{i}.xyz {wd}/{i}_ml')
        os.chdir(f'{wd}/{i}_ml')
        os.system(f'xtb {i}.xyz -c {d[i][0]} --uhf {d[i][1]} --opt --ml_feature>xtb.out')

    dis_gap, dis_pair = spr_intra(smi=list(d.keys())[0],wd=wd,disso=True)
    pull(reactants_temp=list(d.keys()),pull_rec=dis_pair,pull_smi_a=1,wd=wd)

def intermolecular_reaction(target_name_1:str, target_name_2:str, c_1:int, c_2:int, u_1:int, u_2:int):
    wd = os.getcwd()
    d = {}
    d[target_name_1] = [c_1, u_1]
    if target_name_2:
        d[target_name_2] = [c_2, u_2]

    # call xTB
    # os.system('module load xtb/ml-hessian')

    for i in list(d.keys()):
        if os.path.isdir(f'{wd}/{i}_ml'):
            shutil.rmtree(f'{wd}/{i}_ml')
        os.makedirs(f'{wd}/{i}_ml')
        shutil.copy(f'{wd}/{i}.xyz',f'{wd}/{i}_ml')
        os.chdir(f'{wd}/{i}_ml')
        os.system(f'xtb {i}.xyz -c {d[i][0]} --uhf {d[i][1]} --opt --ml_feature>xtb.out')

    # push_dict_2 only for open shell
    push_dict_1, push_dict_2 = spr_inter(list(d.keys()), wd)
    l1 = list(push_dict_1.keys())

    # first two candidates for push_smi_a1 providing HOAO and push_smi_b1 providing LUAO
    push_smi_a1, push_smi_b1 = l1[0].split(",")
    gap_1 = min(sorted(push_dict_1[f'{push_smi_a1},{push_smi_b1}'].keys()))
    gap_2 = sorted(push_dict_1[f'{push_smi_a1},{push_smi_b1}'].keys())[1]
    push_rec_1 = push_dict_1[f'{push_smi_a1},{push_smi_b1}'][gap_1]
    push_rec_2 = push_dict_1[f'{push_smi_a1},{push_smi_b1}'][gap_2]

    # first two candidates for push_smi_a2 providing HOAO and push_smi_b2 providing LUAO
    push_smi_a2, push_smi_b2 = l1[1].split(",")
    gap_3 = min(sorted(push_dict_1[f'{push_smi_a2},{push_smi_b2}'].keys()))
    gap_4 = sorted(push_dict_1[f'{push_smi_a2},{push_smi_b2}'].keys())[1]
    push_rec_3 = push_dict_1[f'{push_smi_a2},{push_smi_b2}'][gap_3]
    push_rec_4 = push_dict_1[f'{push_smi_a2},{push_smi_b2}'][gap_4]

    push(reactants_temp=list(d.keys()), push_rec=push_rec_1, push_smi_a=push_smi_a1, push_smi_b=push_smi_b1, wd=wd)

if __name__ == '__main__':
    if reaction_type == 1:
        call_dissociation()
    if reaction_type == 0:
        intermolecular_reaction()
