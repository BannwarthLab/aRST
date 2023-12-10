import numpy as np
import pandas as pd
from Element import r0


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

        # # Determine the sign of the angle using the dot product between n1 and rc
        if np.dot(n1, rc) < 0.0:
            theta = -theta

        # Convert the angle to degrees
        theta = np.degrees(theta)

        return theta

    def get_ats_distance(self, at1, at2):
        a = self.get_punkt_coord(at1)
        b = self.get_punkt_coord(at2)
        c = b - a
        return np.linalg.norm(c)

    def get_rr_distance(self, at):
        ele = self.coord.iloc[at, 0]
        return r0[ele]


def get_angle(v1, v2):
    r1 = np.round(np.linalg.norm(v1), 1)
    r2 = np.round(np.linalg.norm(v2), 1)
    dot_product = np.round(np.dot(v1, v2), 1)
    try:
        theta = np.arccos(dot_product / (r1 * r2)) * 180 / np.pi
    except RuntimeWarning:
        theta = 90
    return theta


def rerange_P(P, neighbor_atl):
    # determin which main is the reaction direction
    angle_l = []

    if len(neighbor_atl) == 1:
        A = neighbor_atl[0]
        for i in range(3):
            p = np.array(P[:, i])
            angle_l.append(get_angle(p, A))
    else:
        neighbor_atl.append([0, 0, 0])
        A = neighbor_atl[0]
        B = neighbor_atl[1]
        C = neighbor_atl[2]
        AB = np.array(B) - np.array(A)
        AC = np.array(C) - np.array(A)
        normal_vector = np.cross(AB, AC)
        for i in range(3):
            p = np.array(P[:, i])
            dot_product = np.dot(p, normal_vector)
            magnitude_v = np.linalg.norm(p)
            magnitude_normal = np.linalg.norm(normal_vector)
            try:
                theta = np.arccos(dot_product / (magnitude_v * magnitude_normal)) * 180 / np.pi
            except RuntimeWarning:
                theta = 90
            angle_l.append(np.round(theta, 2))

    tmp_l = (abs(np.array(angle_l) - 90)).tolist()
    pr_index = tmp_l.index(max(tmp_l))
    column_order = [0, 1, 2]

    if pr_index < 2:
        for i in range(2 - pr_index):
            column_order.insert(0, column_order.pop())

    P = P[:, column_order]
    return P


def check_orient(P, mv):
    for i in range(3):
        p = np.array(P[:, i])
        theta = get_angle(p, mv)
        if theta < 90:
            P[:, i] = -p
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


def l2r(P, mv):
    # or r2l, it is the same
    theta_l = []
    for i in range(3):
        p = np.array(P[:, i])
        theta = get_angle(p, mv)
        theta_l.append(theta)
    min_index = theta_l.index(min(theta_l))
    tmp = np.array(P[:, min_index])
    P[:, min_index] = -tmp
    return P


def alignment(coord1, coord2, at1, at2, **kwargs):
    # copy coord so that won't affect global variant
    stuc1 = coord1.copy(deep=True)
    stuc2 = coord2.copy(deep=True)

    stuc1 = Coord(stuc1).shift(d=0, at=at1)
    mv1, P1, v1, neighbor_atl1 = Coord(stuc1).get_principal_moment_vec(at1)

    # make sure p point into the direction of the smallest resistance
    P1 = check_orient(P1, mv1)

    # find main reaction direction
    if neighbor_atl1:
        P1 = rerange_P(P1, neighbor_atl1[:3])

    # set P1 left-handed by fliping one main axis to the opposite direction
    # (the one with the least influence on the main orient)
    if not check_rl(P1):
        P1 = l2r(P1, mv1)

    # change coordinante system to main axes system
    stuc1 = Coord(stuc1).change_cs(P1)
    xyz1 = stuc1.copy(deep=True)

    # now do 2. structure
    stuc2 = Coord(stuc2).shift(d=0, at=at2)
    mv2, P2, v2, neighbor_atl2 = Coord(stuc2).get_principal_moment_vec(at2)

    P2 = check_orient(P2, mv2)

    # find main reaction direction
    if neighbor_atl2:
        P2 = rerange_P(P2, neighbor_atl2[:3])

    # set P2 right-handed by fliping one main axis to the opposite direction
    if check_rl(P2):
        P2 = l2r(P2, mv2)

    # test flip all axis? original flip only P2[:, 2]
    P2[:, 0] = -P2[:, 0]
    P2[:, 1] = -P2[:, 1]
    P2[:, 2] = -P2[:, 2]

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

    coord_merge = pd.concat([xyz1, xyz2])
    return coord_merge
