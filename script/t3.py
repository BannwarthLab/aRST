

# stru1_h = [-13.791, -12.7318, -13.8416, -13.8415]
# stru1_l = [-6.7708, -7.9458, 3.3531, 3.3528]
# struc1_c = [0.17363073, -0.27797463, 0.05217271, 0.05217118]
#
# stru2_h = [-14.1855, -13.8043, -13.5864, -15.3651, -13.4637, -13.4637]
# stru2_l = [3.3612, 2.4231, 4.7551, 3.8087, 4.3232, 4.3232]
# struc2_c = [0.0379587, -0.43618072, 0.0569336, 0.28515197, 0.02806822, 0.02806822]
import pandas as pd

# ml1 = pd.read_csv(f'/Users/cy/jtc/cooperation/Buchwald–Hartwig_amination/substrate/xtb/ml/ml_feature.csv')
# ml2 = pd.read_csv(f'/Users/cy/jtc/cooperation/Buchwald–Hartwig_amination/cat/xtb/ml/ml_feature.csv')
# stru1_h = ml1['HOAO_a (eV)'].tolist()
# stru1_l = ml1['LUAO_a (eV)'].tolist()
# struc1_c = pd.read_csv(f'/Users/cy/jtc/cooperation/Buchwald–Hartwig_amination/substrate/xtb/ml/charges',header=None)[0].tolist()
# print(struc1_c[4])
# struc2_c = pd.read_csv(f'/Users/cy/jtc/cooperation/Buchwald–Hartwig_amination/cat/xtb/ml/charges',header=None)[0].tolist()
# print(struc2_c[78])
# stru2_h = ml2['HOAO_a (eV)'].tolist()
# stru2_l = ml2['LUAO_a (eV)'].tolist()

from geom import Coord
from Element import r0, element2
from collections import Counter

import numpy as np

def get_electronic_e(ci, cj, ele1, ele2):

    r = 2 * (r0[ele1] + r0[ele2])
    a = 0
    return a * cj * ci / r

#
# d1 = {}
#
# ml1 = pd.read_csv(f'/Users/cy/jtc/data/formaldehyde/ch3oh/1/1/ml/test/ml_feature_new.csv')
# ml2 = pd.read_csv(f'/Users/cy/jtc/data/formaldehyde/nh3/1/2/ml/test/ml_feature.csv')
# # ml1 = pd.read_csv(f'/Users/cy/jtc/cooperation/Buchwald–Hartwig_amination/substrate/xtb/ml/ml_feature.csv')
# # ml2 = pd.read_csv(f'/Users/cy/jtc/cooperation/Buchwald–Hartwig_amination/cat/xtb/ml/ml_feature.csv')
# ml_df1 = ml1[['atomic partial charges', 'HOAO_a (eV)', 'HOAO_b (eV)', 'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
# ml_df2 = ml2[['atomic partial charges', 'HOAO_a (eV)', 'HOAO_b (eV)', 'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom']]
# struc1_c = ml_df1['atomic partial charges'].tolist()
# struc2_c = ml_df2['atomic partial charges'].tolist()
# ele_l1 = ml_df1['Atom'].tolist()
# ele_l2 = ml_df2['Atom'].tolist()
# h1 = ml_df1['HOAO_a (eV)'].tolist()
# h2 = ml_df2['HOAO_a (eV)'].tolist()
# l1 = ml_df1['LUAO_a (eV)'].tolist()
# l2 = ml_df2['LUAO_a (eV)'].tolist()
# pair_rec = []
# for i in range(len(h1)):
#     for j in range(len(h2)):
#         ele1 = element2[ele_l1[i]]
#         ele2 = element2[ele_l2[j]]
#         ci = struc1_c[i]
#         cj = struc2_c[j]
#         e = get_electronic_e(ci, cj, ele1, ele2)
#         key = np.round(l2[j] - h1[i] + e, 3)
#         d1[key] = f'{i + 1},{j + 1},h,gap = {np.round(l2[j] - h1[i], 3)},e = {np.round(e, 3)}'
#         pair_rec.append([i + 1, j + 1]) if not any(Counter([i + 1, j + 1]) == Counter(k) for k in pair_rec) else None
#
#         key = np.round(l1[i] - h2[j] + e, 3)
#         d1[key] = f'{i + 1},{j + 1},l,gap = {np.round(l1[i] - h2[j], 3)},e = {np.round(e, 3)}'
#         pair_rec.append([i + 1, j + 1]) if not any(Counter([i + 1, j + 1]) == Counter(k) for k in pair_rec) else None
# d1 = {key: value for key, value in d1.items() if value is not None}
#
# print(len(d1.keys()))
# for i in list(sorted(d1.keys()))[:30]:
#     print(i, d1[i])

# ml1 = pd.read_csv(f'/Users/cy/jtc/data/formaldehyde/ch3oh/1/1/ml/test/ml_feature_new.csv')
# ml2 = pd.read_csv(f'/Users/cy/jtc/data/formaldehyde/nh3/1/2/ml/test/ml_feature.csv')
# HOAO_a = pd.DataFrame()
# for i in[ml1,ml2]:
#     HOAO_a = pd.concat([HOAO_a, i['HOAO_a (eV)']], axis=1)
# print(HOAO_a)

# ml = pd.read_csv(f'/Users/cy/jtc/data/formaldehyde/ch3oh/struc/0001/ml/ml_feature.csv', sep=',')
#
# for i in ml.columns:
#     print(ml[i])


def get_punkt_coord(coord,at):
    a = np.array([coord.iloc[at - 1, 1],
                  coord.iloc[at - 1, 2],
                  coord.iloc[at - 1, 3]])

    return a
def get_dihedral_angle(coord, at1, at2, at3, at4):
    a = get_punkt_coord(coord, at1)
    b = get_punkt_coord(coord, at2)
    c = get_punkt_coord(coord, at3)
    d = get_punkt_coord(coord, at4)

    ra = b - a
    rb = c - b
    rc = d - c

    n1 = np.cross(ra, rb)
    n2 = np.cross(rb, rc)

    n1_normalized = n1 / np.linalg.norm(n1)
    n2_normalized = n2 / np.linalg.norm(n2)

    # Calculate the dihedral angle in radians
    # sign = np.sign(np.dot(np.cross(ra, rb), rc))
    theta = np.arccos(np.dot(n1_normalized, n2_normalized))

    # # Determine the sign of the angle using the dot product between n1 and rc
    if np.dot(n1, rc) < 0.0:
        theta = -theta

    # Convert the angle to degrees
    theta = np.degrees(theta)

    return sign

# at1 = 25
# at2 = 13
# at3 = 3
# at4 = 24
# coord = pd.read_table('coord.xyz', header=None, skiprows=2, sep='\s+')
# print(get_dihedral_angle(coord,at1,at2,at3,at4))
# print(get_dihedral_angle(coord,at1,at3,at2,at4))
reactants_intraa, products_intraa = [[1],[2]],[[1],[2]]

mllevel = "dft_pbe"
if "gfn" in mllevel or mllevel is None:
#
# if not "dft" in mllevel or not "DFT" in mllevel or mllevel is None:
    print("a")
else:
    print("b")

