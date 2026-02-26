import os
import numpy as np
from aRST.script.element import r0, element
import warnings
import math

warnings.filterwarnings("ignore")


class Coord:
    def __init__(self, coordinate, ele_l, centerats=None):
        self.coordinate = coordinate
        self.ele_l = ele_l
        self.rcov_l = np.array([float(r0[e]) for e in self.ele_l], dtype=np.float64)
        self.rcov_l_scaled = np.array([(4 / 3) * _ for _ in self.rcov_l], dtype=np.float64)
        self.nat = len(self.ele_l)
        self.nat_l = np.array([float(element[e]) for e in self.ele_l], dtype=np.int16)
        self.centerats = centerats

    def get_principal_moment_vec(self, centerats_ext=None):
        t1, t2, t3, t4, t5, t6 = 0, 0, 0, 0, 0, 0
        # matrix for moments of inertia is of form, with CN in mass as weight
        # *           |   y**2+z**2                         |
        # *           |    -y*x       z**2+x**2             | -i =0
        # *           |    -z*x        -z*y       x**2+y**2 |
        xv, yv, zv = 0, 0, 0
        neighbor_atl = []

        if centerats_ext is None:
            coords = [self.coordinate[i] for i in self.centerats]
            ax, ay, az = [sum(values) / len(values) for values in zip(*coords)]
            rcov_c_scaled = sum(self.rcov_l_scaled[i] for i in self.centerats) / len(self.centerats)
            centerats = self.centerats
        else:
            coords = [self.coordinate[i] for i in centerats_ext]
            ax, ay, az = [sum(values) / len(values) for values in zip(*coords)]
            rcov_c_scaled = sum(self.rcov_l_scaled[i] for i in centerats_ext) / len(centerats_ext)
            centerats = centerats_ext

        for i in range(len(self.ele_l)):
            if i in centerats:
                continue
            bx, by, bz = self.coordinate[i]
            x = bx - ax
            y = by - ay
            z = bz - az
            r2 = x * x + y * y + z * z
            r = math.sqrt(r2)
            rcov_i_scaled = self.rcov_l_scaled[i]
            rr = (rcov_i_scaled + rcov_c_scaled) / r

            # neighbors
            if rr > 0.8:
                neighbor_atl.append([bx, by, bz])
            cn = 1.0 / (1.0 + math.exp(-16.0 * (rr - 1.0)))
            atmass = 0.5 * cn / r2

            # inertia tensor components
            t1 += atmass * (y * y + z * z)
            t2 += atmass * (-y * x)
            t3 += atmass * (z * z + x * x)
            t4 += atmass * (-z * x)
            t5 += atmass * (-z * y)
            t6 += atmass * (x * x + y * y)

            # add up molecular orientation
            w = self.nat_l[i]
            xv += w * x
            yv += w * y
            zv += w * z
        mv = np.array([xv, yv, zv], dtype=np.float64)

        arep = np.array([
            [t1, t2, t4],
            [t2, t3, t5],
            [t4, t5, t6]
        ], dtype=np.float64)

        v, P = np.linalg.eigh(arep)
        #     print(v.round(3))
        #     print(P.round(3))
        return mv, P, v, neighbor_atl

    def shift(self, d, centerats_ext=None):
        # shift the molecule with aktiv as center to (0,0,d)
        if centerats_ext is None:
            coords = [self.coordinate[i] for i in self.centerats]
            shift_vec = np.array([sum(values) / len(values) for values in zip(*coords)])
        else:
            coords = [self.coordinate[i] for i in centerats_ext]
            shift_vec = np.array([sum(values) / len(values) for values in zip(*coords)])

        shift_vec[2] -= d
        self.coordinate -= shift_vec
        return self

    def change_cs(self, P):
        self.coordinate = self.coordinate @ P
        return self

    def write_xyz(self, name):
        with open(f"{name}.xyz", "w") as f:
            f.write(f"{self.nat}\n\n")
            for at, (x, y, z) in zip(self.ele_l, self.coordinate):
                f.write(f"{at:2s} {x:18.14f} {y:18.14f} {z:18.14f}\n")
            f.write("\n")

    def write_p(self, P, name, dst, length=2, porder=None):
        # display p vector in chimera, defalut length =2
        coords = [self.coordinate[i] for i in self.centerats]
        x, y, z = [sum(values) / len(values) for values in zip(*coords)]
        P_scalred = np.copy(P)
        for i in range(P_scalred.shape[1]):
            tmp = np.linalg.norm(P[:, i])
            P_scalred[:, i] = P_scalred[:, i] * length / tmp

        targetfile = os.path.join(dst, f'{name}.bild')
        color0 = ["green", "blue", "red"]
        if porder is not None:
            color = [color0[i] for i in porder]
        else:
            color = color0
        with open(targetfile, 'w+') as f2:
            # Coordinate system
            f2.write('.color black\n')
            f2.write(f'.arrow 0 0 0 0 0 6 0.05 0.150000\n')
            f2.write(f'.arrow 0 0 0 6 0 0 0.05 0.150000\n')
            f2.write(f'.arrow 0 0 0 0 6 0 0.05 0.150000\n')

            f2.write(f'.color {color[2]}\n')
            f2.write(
                f'.arrow {x} {y} {z} {x + P_scalred[:, 2][0]} {y + P_scalred[:, 2][1]} {z + P_scalred[:, 2][2]} 0.05 0.150000\n')

            f2.write(f'.color {color[0]}\n')
            f2.write(
                f'.arrow {x} {y} {z} {x + P_scalred[:, 0][0]} {y + P_scalred[:, 0][1]} {z + P_scalred[:, 0][2]} 0.05 0.150000\n')

            f2.write(f'.color {color[1]}\n')
            f2.write(
                f'.arrow {x} {y} {z} {x + P_scalred[:, 1][0]} {y + P_scalred[:, 1][1]} {z + P_scalred[:, 1][2]} 0.05 0.150000\n')
        return

    def get_dihedral_angle(self, at1, at2, at3, at4):
        # struc1: at1, at3 and at1 = 0,0,0
        # struc2: at2,at4 and at2 = 0,0,z
        a = self.coordinate[at1]
        b = self.coordinate[at2]
        c = self.coordinate[at3]
        d = self.coordinate[at4]

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

    def choose_combination(self, at1, at2, at3, at4):
        # atl =  [[at1,at2],[at3,at4]]
        # at1 = 0,0,0; at2 = 0,0,z
        d12 = self.get_ats_distance(at1, at2)
        d34 = self.get_ats_distance(at3, at4)
        d14 = self.get_ats_distance(at1, at4)
        d32 = self.get_ats_distance(at3, at2)
        straight = d12 + d34
        swapped = d14 + d32
        if straight <= swapped:
            return "straight"
        else:
            return "swapped"

    def get_ats_distance(self, at1, at2):
        ax, ay, az = self.coordinate[at1]
        bx, by, bz = self.coordinate[at2]

        x = bx - ax
        y = by - ay
        z = bz - az

        r2 = x * x + y * y + z * z
        r = math.sqrt(r2)
        return r


def get_angle(v1, v2):  # 0-180
    r1 = np.round(np.linalg.norm(v1), 3)
    r2 = np.round(np.linalg.norm(v2), 3)
    dot_product = np.round(np.dot(v1, v2), 3)
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
    column_order = [0, 1, 2]
    if v[1] == v[2] and v[0] != 0:  # prolate
        column_order = [0, 1, 2]
        pr_index = v.index(min(v))
        for i in range(2 - pr_index):
            column_order.insert(0, column_order.pop())
        P = P[:, column_order]
    else:
        P = P
    return P, column_order


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


def alignment(coordinate1, ele_l1, coordinate2, ele_l2, at_l, extra_setting=None, writingP=False):
    # copy coord so that won't affect global variant
    struc1 = coordinate1.copy()
    struc2 = coordinate2.copy()

    centerats_1 = [int(_[0]) for _ in at_l]
    struc1 = Coord(struc1, ele_l1, centerats=centerats_1).shift(d=0)
    mv1, P1, v1, neighbor_atl1 = struc1.get_principal_moment_vec()

    # make sure p point into the direction of the smallest resistance and keep right-handed
    P1 = check_orient(P1, mv1)
    if writingP:
        struc1.write_xyz("coord1")
        struc1.write_p(P1, "p1_primitive", os.getcwd(), 2)

    # find main reaction direction
    P1, new_porder = check_range(P1, v1)

    # change coordinate system to main axes system
    struc1.change_cs(P1)
    if writingP:
        struc1.write_p(P1, "p1_oriented", os.getcwd(), 2, new_porder)

    # now do 2. structure
    centerats_2 = [int(_[1]) for _ in at_l]
    struc2 = Coord(struc2, ele_l2,centerats=centerats_2).shift(d=0)
    mv2, P2, v2, neighbor_atl2 = struc2.get_principal_moment_vec()

    P2 = check_orient(P2, mv2)
    if writingP:
        struc2.write_xyz("coord2")
        struc2.write_p(P2, "p2_primitive", os.getcwd(), 2)

    # find main reaction direction
    P2, new_porder = check_range(P2, v2)

    P2[:, 2] = -P2[:, 2]

    P2[:, 0] = -P2[:, 0]
    P2[:, 1] = -P2[:, 1]

    maxovlap = extra_setting.get("maxovlap", False)
    if maxovlap:
        P2[:, 0] = -P2[:, 0]
        P2[:, 1] = -P2[:, 1]

    struc2.change_cs(P2)
    if writingP:
        struc2.write_p(P2, "p2_oriented", os.getcwd(), 2, new_porder)

    # shift struc2 to (0,0,d)
    cov_r_at12 = r0[struc1.ele_l[centerats_1[0]]] + r0[struc2.ele_l[centerats_2[0]]]
    relativ_r = extra_setting.get('aligning_distance', 5)  # 5 times of rr by default
    d = np.round(cov_r_at12 * relativ_r, 3)
    struc2 = struc2.shift(d=d)

    # merge two structures
    coordinate_merge = np.vstack((struc1.coordinate, struc2.coordinate))
    ele_merge = np.concatenate([struc1.ele_l, struc2.ele_l])

    return coordinate_merge, ele_merge
