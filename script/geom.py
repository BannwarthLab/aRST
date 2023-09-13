import numpy as np
import pandas as pd

r0 = {'H': 0.32, 'D': 0.32, 'He': 0.46, 'Li': 1.20, 'Be': 0.94, 'B': 0.77, 'C': 0.75, 'N': 0.71, 'O': 0.63,
      'F': 0.64, 'Ne': 0.67, 'Na': 1.40, 'Mg': 1.25, 'Al': 1.13, 'Si': 1.04, 'P': 1.10, 'S': 1.02,
      'CL': 0.99, 'Ar': 0.96, 'K': 1.76, 'Ca': 1.54, 'Sc': 1.48, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22,
      'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11, 'Ni': 1.10, 'Cu': 1.12, 'Zn': 1.18, 'Ga': 1.24, 'Ge': 1.21,
      'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17, 'Rb': 2.10, 'Sr': 1.85, 'Y': 1.63, 'Zr': 1.54,
      'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25, 'Pd': 1.20, 'Ag': 1.28, 'Cd': 1.36,
      'In': 1.42, 'Sn': 1.40, 'Sb': 1.40, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31, 'Cs': 2.32, 'Ba': 1.96,
      'La': 1.80, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74, 'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68, 'Gd': 1.69,
      'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.70, 'Lu': 1.62, 'Hf': 1.52,
      'Ta': 1.46, 'W': 0.95 * 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24, 'Hg': 1.33,
      'Tl': 1.44, 'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01,
      'Ac': 1.86, 'Th': 1.75, 'Pa': 1.69, 'U': 1.70, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66,
      'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Rf': 1.57,
      'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34, 'Mt': 1.29, 'Ds': 1.21, 'Rg': 1.21, 'Cn': 1.22,
      'Nh': 1.36, 'Fl': 1.43, 'Mc': 1.62, 'Lv': 1.75, 'Ts': 1.65, 'Og': 1.57, 'Cl': 0.99}


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
        r1 = np.round(np.linalg.norm(np.array(v1)), 1)
        r2 = np.round(np.linalg.norm(np.array(v2)), 1)
        dot_p = np.round(np.dot(v1, v2), 1)
        c = np.round(dot_p / (r1 * r2), 1)
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
            # a = [x, y, z]
            a_new = [x, -y, z]
            self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3] = a_new[0], a_new[1], a_new[2]
        return self.coord

    def get_dihral_angle(self, at1, at2, at3, at4):
        a = Coord.get_punkt_coord(self, at1)
        b = Coord.get_punkt_coord(self, at2)
        c = Coord.get_punkt_coord(self, at3)
        d = Coord.get_punkt_coord(self, at4)
        print(a, b, c, d)
        ra = b - a
        rb = c - b
        rc = d - c
        n1 = np.cross(ra, rb)
        n2 = np.cross(rb, rc)
        # arccos（n1·n2）/（|n1|*|n2|)
        theta = np.rad2deg(np.arccos(np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))))
        return theta


def alignment(coord1, coord2, at1, at2, **kwargs):
    C1 = Coord(coord1)
    coord1 = C1.shift(at=at1, d=0)
    name1 = coord1.iloc[at1 - 1, 0]
    rhanded = C1.test_if_right_handed_system(at=at1)
    orient = C1.geom_orientation(at=at1)
    P1 = C1.get_principal_moment_vec(at=at1)[1]
    if not rhanded:
        coord1 = C1.l2r()  # flip y to -y
    if not orient:
        coord1 = C1.rot_x()  # rot around x-axis 180
    coord1 = Coord(coord1).change_cs(P1).copy(deep=True)  # for the bug if cood1=cood2, the class will be inherited

    C2 = Coord(coord2)
    coord2 = C2.shift(at=at2, d=0)
    name2 = coord2.iloc[at2 - 1, 0]
    rhanded = C2.test_if_right_handed_system(at=at2)
    orient = C2.geom_orientation(at=at2)
    P2 = C2.get_principal_moment_vec(at=at2)[1]
    if not rhanded:
        coord2 = C2.l2r()  # flip y to -y
    if not orient:  # Pz point into together
        coord2 = C2.rot_x()  # rot around x-axis 180
    coord2 = Coord(coord2).change_cs(P2)

    rr = (r0[name1] + r0[name2])
    relativ_r = kwargs.get('aligning_distance', 5)  # 5 times of rr by default
    d_pre = kwargs.get('aligning_direction', 1)  # 1 by default, -1: opposite alignment
    d = d_pre * np.round(rr * relativ_r, 2)
    coord2 = Coord(coord2).shift(at=at2, d=d).copy(deep=True)

    coord_merge = pd.concat([coord1, coord2])
    return coord_merge, rr
