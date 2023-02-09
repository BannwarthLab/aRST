import os
import shutil

import numpy as np
import pandas as pd

"""Please place this script and your .xyz file in the same folder

please load module : xtb/ml-hessian !!!!!


if you want to predict the intermolecular reaction, fill in the target_name_1/2 
e.g. target_name_1 = str("ch3br"), target_name_2 = str("cl").

in the last line you can also change different reaction sites candidates
"""

target_name_1 = str("a")
c_1 = 0
u_1 = 0  # !!! please change here if charge and uhf are not 0!
target_name_2 = str("b")
c_2 = 0
u_2 = 0  # !!! please change here if charge and uhf are not 0!

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


#######################Functions#########################

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
        d = {}
        d_1 = {}
        d_2 = {}
        for i in range(len(h)):
            if i % 2 == 0:
                start = 1
            else:
                start = 0
            for j in range(*{'start': start, 'stop': len(l), 'step': 2}.values()):
                if int(i / 2) != int(j / 2):
                    h_i = h[i]
                    l_j = l[j]
                    # print(h_i)
                    # print(l_j)
                    d[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = {}
                    d_2[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = {}
                    for one in range(len(h_i)):
                        for two in range(len(l_j)):
                            try:
                                h_i[one] = float(h_i[one])
                            except ValueError:
                                h_i[one] = float('inf')
                            try:
                                l_j[two] = float(l_j[two])
                            except ValueError:
                                l_j[two] = float('inf')

                            key = np.round(abs(float(h_i[one]) - float(l_j[two])), 3)
                            if key in d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].keys():
                                d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key].append([one + 1, two + 1])
                            else:
                                d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key] = [[one + 1, two + 1]]

                            for item_2 in list(d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].keys()):
                                if 0 < abs(key - item_2) < 0.01:
                                    temp = list(d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][item_2])
                                    temp[0].append('*')
                                    d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key] = \
                                        d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][
                                            key] + temp
                                    del d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][item_2]
                    # print(d[f'{int(i/2)+1},{int(j/2)+1}'])
                    # print(int(i/2)+1, int(j/2)+1)
                    # Discard excess data and organize it
                    if f'{int(i / 2) + 1},{int(j / 2) + 1}' in d_1.keys():
                        d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] + sorted(
                            d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].items())
                    else:
                        d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = sorted(
                            d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].items())
                    # print(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'])
                    if len(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}']) > 1:
                        if d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][0][0] < \
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][1][
                                    0]:
                            d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'].remove(
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][1])
                        else:
                            d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'].remove(
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][0])
                    for u in list(dict(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}']).keys()):
                        l_1 = []
                        l_2 = []
                        for item in list(dict(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'])[u]):
                            if (item[0] not in l_2) and (item[1] not in l_2):
                                l_2.append(item[0])
                                l_2.append(item[1])
                                l_1.append(item)
                                d_2[f'{int(i / 2) + 1},{int(j / 2) + 1}'][u] = l_1
                                # dict(d_1[f'{i + 1},{j + 1}'])[u] = l_1
        for a in d_2.keys():
            for b in list(d_2[a].keys())[2:]:
                del d_2[a][b]
        # print("push_dict: alpha beta: ")
        # print(d_2)

        d = {}
        d_1 = {}
        d_3 = {}
        for i in range(len(h)):
            if i % 2 != 0:
                start = 1
            else:
                start = 0
            for j in range(*{'start': start, 'stop': len(l), 'step': 2}.values()):
                if int(i / 2) != int(j / 2):
                    h_i = h[i]
                    l_j = l[j]
                    # print(h_i)
                    # print(l_j)
                    d[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = {}
                    d_3[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = {}
                    for one in range(len(h_i)):
                        for two in range(len(l_j)):
                            try:
                                h_i[one] = float(h_i[one])
                            except ValueError:
                                h_i[one] = float('inf')
                            try:
                                l_j[two] = float(l_j[two])
                            except ValueError:
                                l_j[two] = float('inf')
                            key = np.round(abs(float(h_i[one]) - float(l_j[two])), 3)
                            if key in d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].keys():
                                d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key].append([one + 1, two + 1])
                            else:
                                d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key] = [[one + 1, two + 1]]

                            for item_2 in list(d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].keys()):
                                if 0 < abs(key - item_2) < 0.01:
                                    temp = list(d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][item_2])
                                    temp[0].append('*')
                                    d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][key] = \
                                        d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][
                                            key] + temp
                                    del d[f'{int(i / 2) + 1},{int(j / 2) + 1}'][item_2]

                    if f'{int(i / 2) + 1},{int(j / 2) + 1}' in d_1.keys():
                        d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] + sorted(
                            d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].items())
                    else:
                        d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'] = sorted(
                            d[f'{int(i / 2) + 1},{int(j / 2) + 1}'].items())
                    # print(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'])
                    if len(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}']) > 1:
                        if d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][0][0] < \
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][1][
                                    0]:
                            d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'].remove(
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][1])
                        else:
                            d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'].remove(
                                d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'][0])
                    for u in list(dict(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}']).keys()):
                        l_1 = []
                        l_2 = []
                        for item in list(dict(d_1[f'{int(i / 2) + 1},{int(j / 2) + 1}'])[u]):
                            if (item[0] not in l_2) and (item[1] not in l_2):
                                l_2.append(item[0])
                                l_2.append(item[1])
                                l_1.append(item)
                                d_3[f'{int(i / 2) + 1},{int(j / 2) + 1}'][u] = l_1
                                # dict(d_1[f'{i + 1},{j + 1}'])[u] = l_1
        for a in d_3.keys():
            for b in list(d_3[a].keys())[2:]:
                del d_3[a][b]
        # print("push_dict: alpha alpha: ")
        # print(d_3)
        print(d_2, d_3)
        return d_2, d_3
    else:
        d_2 = {}
        d_3 = {}
        return d_2, d_3


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


def aligment(smi_a, smi_b, rec, d, wd):
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
        f2.write(f'.arrow 0 0 0 {P1[2][0]*4} {P1[2][1]*4} {P1[2][2]*4} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 0 {P1[0][0]} {P1[0][1]} {P1[0][2]} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 0 {P1[1][0]*2} {P1[1][1]*2} {P1[1][2]*2} 0.05 0.150000\n')


        f2.write('.color green\n')
        f2.write(f'.arrow 0 0 0 10 0 0 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[2][0]*4} {P2[2][1]*4} {(P2[2][2]*4 + d)} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[0][0]} {P2[0][1]} {P2[0][2] + d} 0.05 0.150000\n')
        f2.write(f'.arrow 0 0 {d} {P2[1][0]*2} {P2[1][1]*2} {(P2[1][2]*2 + d)} 0.05 0.150000\n')

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
    atomnnr_a = aligment(smi_a=push_smi_a, smi_b=push_smi_b, rec=push_rec[0], d=align_distance, wd=wd)
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


#########################################################

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
    os.system(f'cp {wd}/{i}.xyz {wd}/{i}_ml')
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

# you can also check other candidates:
#e.g. :{'1,2': {10.325000000000001: [[2, 2]], 13.029: [[1, 1]]}, '2,1': {10.653: [[1, 2]], 12.012: [[2, 1]]}}
# push_dict_2: first molecule provides HOAO and second provides LUAO, [1, 1]
# push_dict_3: first molecule provides LUAO and second provides HOAO , [1, 2]
# push_dict_4: first molecule provides LUAO and second provides HOAO , [2, 1]
# then please change last line to
# push(reactants_temp=list(d.keys()), push_rec=push_rec_2, push_smi_a=push_smi_a1, push_smi_b=push_smi_b1, wd=wd)
# push(reactants_temp=list(d.keys()), push_rec=push_rec_3, push_smi_a=push_smi_a2, push_smi_b=push_smi_b2, wd=wd)
# push(reactants_temp=list(d.keys()), push_rec=push_rec_4, push_smi_a=push_smi_a2, push_smi_b=push_smi_b2, wd=wd)

