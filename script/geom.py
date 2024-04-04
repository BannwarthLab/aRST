import numpy as np
import pandas as pd
from Element import r0,element


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
        neighbor_atl = []
        for i in range(self.coord.shape[0]):
            if i != at:
                x = self.coord.iloc[i, 1] - self.coord.iloc[at, 1]
                y = self.coord.iloc[i, 2] - self.coord.iloc[at, 2]
                z = self.coord.iloc[i, 3] - self.coord.iloc[at, 3]
                r = np.linalg.norm(np.array([x, y, z]))
                rcov_c = float(r0[f'{self.coord.iloc[at, 0]}'])
                rcov_i = float(r0[f'{self.coord.iloc[i, 0]}'])
                rr = (rcov_i + rcov_c) / r
                # add connected at's xyz
                if rr > 0.8:
                    neighbor_atl.append([self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3]])
                # add up I
                cn = 1 / (1 + np.exp(- 16 * rr - 1))
                atmass = 0.5 * cn / (r ** 2)
                t1 = t1 + atmass * (y ** 2 + z ** 2)
                t2 = t2 + atmass * (-y * x)
                t3 = t3 + atmass * (z ** 2 + x ** 2)
                t4 = t4 + atmass * (-z * x)
                t5 = t5 + atmass * (-z * y)
                t6 = t6 + atmass * (x ** 2 + y ** 2)
                # add up molecular orientation
                xv = xv + element[f'{self.coord.iloc[i, 0]}']* x
                yv = yv + element[f'{self.coord.iloc[i, 0]}'] * y
                zv = zv + element[f'{self.coord.iloc[i, 0]}'] * z
        # print(cn)
        mv = np.array([xv, yv, zv])
        arep = np.zeros((3, 3))
        arep[0][0], arep[1][0], arep[1][1], arep[2][0], arep[2][1], arep[2][2] = t1, t2, t3, t4, t5, t6
        arep = arep + arep.T - np.diag(arep.diagonal())

        v, P = np.linalg.eigh(arep)
        #     print(v.round(3))
        #     print(P.round(3))
        return mv, P, v, neighbor_atl

    def get_punkt_coord(self, at):
        a = np.array([self.coord.iloc[at, 1],
                      self.coord.iloc[at, 2],
                      self.coord.iloc[at, 3]])
        return a

    def get_nat(self):
        return self.coord.shape[0]

    def shift(self, at, d):
        # shift the molecule with aktiv as center to (0,0,d)
        x = self.coord.iloc[at, 1]
        y = self.coord.iloc[at, 2]
        z = self.coord.iloc[at, 3] - d
        for i in range(self.coord.shape[0]):
            self.coord.iloc[i, 1] = self.coord.iloc[i, 1] - x
            self.coord.iloc[i, 2] = self.coord.iloc[i, 2] - y
            self.coord.iloc[i, 3] = self.coord.iloc[i, 3] - z
        return self.coord

    def change_cs(self, P):
        for i in range(self.coord.shape[0]):
            x, y, z = self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3]
            a = [x, y, z]
            a_new = np.dot(P.T, a)
            self.coord.iloc[i, 1], self.coord.iloc[i, 2], self.coord.iloc[i, 3] = a_new[0], a_new[1], a_new[2]
        return self.coord

    def write_xyz(self, name):
        xyz = self.coord.round(5)
        nat = xyz.shape[0]
        with open(f'{name}.xyz', 'w') as f:
            f.write(f'{nat}\n')
            f.write('\n')
            for i in range(nat):
                at, x, y, z = xyz.iloc[i, 0], xyz.iloc[i, 1], xyz.iloc[i, 2], xyz.iloc[i, 3]
                f.write('%2s %10.5f %10.5f %10.5f\n' % (at, x, y, z))
            f.write('\n')

    def write_p(self, at, P, v, name):
        x, y, z = 0, 0, 0
        try:
            a = Coord.get_punkt_coord(self, at)
            x, y, z = a[0], a[1], a[2]
        except:
            pass
        if P is not None and v is not None:
            scale = 0.5
            with open(f'{name}.bild', 'w+') as f2:
                f2.write('.color black\n')
                f2.write(f'.arrow 0 0 0 0 0 6 0.05 0.150000\n')
                f2.write(f'.arrow 0 0 0 6 0 0 0.05 0.150000\n')
                f2.write(f'.arrow 0 0 0 0 6 0 0.05 0.150000\n')

                f2.write('.color red\n')
                f2.write(
                    f'.arrow {x} {y} {z} {x + P[::, 2][0] / scale} {y + P[::, 2][1] / scale} {z + P[::, 2][2] / scale} 0.05 0.150000\n')
                f2.write(
                    f'.arrow {x} {y} {z} {x - P[::, 2][0] / scale} {y - P[::, 2][1] / scale} {z - P[::, 2][2] / scale} 0.05 0.150000\n')

                f2.write('.color green\n')
                f2.write(
                    f'.arrow {x} {y} {z} {x + P[::, 0][0] / scale} {y + P[::, 0][1] / scale} {z + P[::, 0][2] / scale} 0.05 0.150000\n')
                f2.write(
                    f'.arrow {x} {y} {z} {x - P[::, 0][0] / scale} {y - P[::, 0][1] / scale} {z - P[::, 0][2] / scale} 0.05 0.150000\n')

                f2.write('.color blue\n')
                f2.write(
                    f'.arrow {x} {y} {z} {x + P[::, 1][0] / scale} {y + P[::, 1][1] / scale} {z + P[::, 1][2] / scale} 0.05 0.150000\n')
                f2.write(
                    f'.arrow {x} {y} {z} {x - P[::, 1][0] / scale} {y - P[::, 1][1] / scale} {z - P[::, 1][2] / scale} 0.05 0.150000\n')
        else:
            with open(f'{name}.bild', 'w+') as f2:
                f2.write('.color black\n')
                f2.write(f'.arrow 0 0 0 0 0 6 0.05 0.150000\n')
                f2.write(f'.arrow 0 0 0 6 0 0 0.05 0.150000\n')
                f2.write(f'.arrow 0 0 0 0 6 0 0.05 0.150000\n')
        return

    def get_dihedral_angle(self, at1, at2, at3, at4):
        # struc1: at1, at3 and at1 = 0,0,0
        # struc2: at2,at4 and at2 = 0,0,z
        a = self.get_punkt_coord(at1)
        b = self.get_punkt_coord(at2)
        c = self.get_punkt_coord(at3)
        d = self.get_punkt_coord(at4)

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

        #         # # Determine the sign of the angle using the dot product between n1 and rc
        if np.dot(n1, rc) < 0.0:
            theta = -theta

        # Convert the angle to degrees
        theta = np.degrees(theta)

        return theta

    def test_if_cross(self, at1, at2, at3, at4):
        # atl =  [[at1,at2],[at3,at4]]
        # at1 = 0,0,0; at2 = 0,0,z
        d1 = self.get_dihedral_angle(at3, at1, at2, at4)
        d2 = self.get_dihedral_angle(at3, at1, at4, at2)

        if abs(d1) >= abs(d2):
            return True
        else:
            return False

    def get_ats_distance(self, at1, at2):
        a = self.get_punkt_coord(at1)
        b = self.get_punkt_coord(at2)
        c = b - a
        return np.linalg.norm(c)

    def get_rr_distance(self, at):
        ele = self.coord.iloc[at, 0]
        return r0[ele]


def get_angle(v1, v2): # 0-180
    r1 = np.round(np.linalg.norm(v1), 1)
    r2 = np.round(np.linalg.norm(v2), 1)
    dot_product = np.round(np.dot(v1, v2), 1)
    try:
        theta = np.arccos(dot_product / (r1 * r2)) * 180 / np.pi
    except RuntimeWarning:
        theta = 90
    return theta


def check_orient(P, mv, rh_cs=True):
    # default to optput right-handed P
    tmp, tmp2 = [], []
    current_type = check_rl(P)  # true = left handed

    for i in range(3):
        p = np.array(P[:, i])
        theta = get_angle(p, mv)
        tmp2.append(abs(abs(theta) - 90))

        if theta < 90:
            tmp.append(i)
            P[:, i] = -p

    if len(tmp) == 1 or len(tmp) == 3:
        if (rh_cs and current_type) or (not rh_cs and not current_type):
            i = tmp2.index(min(tmp2))
            P[:, i] = -P[:, i]
        else:
            return P

    return P


def check_range(P, v):
    #
    v = v.round(2)
    v = [i for i in v]
    # a=b=c spherical pr=pc
    # a<b=c prolate pr=pa !
    # a=b<c oblate pr=pc
    # a=0 b=c linear pr=pc
    # a+b = c oblate, plane pr=pc

    if v[1] == v[2] and v[0] != 0:  # prolate
        column_order = [0, 1, 2]
        pr_index = v.index(min(v))
        for i in range(2 - pr_index):
            column_order.insert(0, column_order.pop())
        P = P[:, column_order]
    else:
        P = P
    return P


def check_rl(P):
    # retrun true if P is left-handed
    a, b, c = P[::, 0], P[::, 1], P[::, 2]
    cross_p = [np.round(j, 3) for j in np.cross(b, c)]
    dot_p = np.dot(a, cross_p)
    if dot_p > 0:
        return True
    else:
        return False


def l2r(P):
    column_order = [1, 0, 2]
    P = P[:, column_order]
    return P


def alignment(coord1, coord2, at1, at2, **kwargs):
    # copy coord so that won't affect global variant
    stuc1 = coord1.copy(deep=True)
    stuc2 = coord2.copy(deep=True)

    stuc1 = Coord(stuc1).shift(d=0, at=at1)
    mv1, P1, v1, neighbor_atl1 = Coord(stuc1).get_principal_moment_vec(at1)

    # make sure p point into the direction of the smallest resistance and keep right-handed
    P1 = check_orient(P1, mv1)


    # find main reaction direction
    P1 = check_range(P1,v1)


    # change coordinante system to main axes system
    stuc1 = Coord(stuc1).change_cs(P1)
    xyz1 = stuc1.copy(deep=True)

    # now do 2. structure
    stuc2 = Coord(stuc2).shift(d=0, at=at2)
    mv2, P2, v2, neighbor_atl2 = Coord(stuc2).get_principal_moment_vec(at2)

    P2 = check_orient(P2, mv2)

    # find main reaction direction
    P2 = check_range(P2,v2)

    P2[:, 2] = -P2[:, 2]

    P2[:, 0] = -P2[:, 0]
    P2[:, 1] = -P2[:, 1]

    ovlap_mode = kwargs.get('aligning_ovlap', False)
    if ovlap_mode:
        P2[:, 0] = -P2[:, 0]
        P2[:, 1] = -P2[:, 1]

    stuc2 = Coord(stuc2).change_cs(P2)

    # merge two structures
    name1 = xyz1.iloc[at1, 0]
    name2 = stuc2.iloc[at2, 0]
    rr = (r0[name1] + r0[name2])
    relativ_r = kwargs.get('aligning_distance', 5)  # 5 times of rr by default
    d_pre = kwargs.get('aligning_direction', 1)  # 1 by default, -1: opposite alignment
    d = d_pre * np.round(rr * relativ_r, 2)

    stuc2 = Coord(stuc2).shift(at2, d)
    xyz2 = stuc2.copy(deep=True)  # for the bug if cood1=cood2, the class will be inherited

    coord_merge = pd.concat([xyz1, xyz2], ignore_index=True)
    return coord_merge
