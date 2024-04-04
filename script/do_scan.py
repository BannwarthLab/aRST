import os

import numpy as np
import pandas as pd
from toolbox import create_folder, write_xyz, get_nelectron
from spr import SPR_HL, build_bo_d
from molbar.barcode import idealize_structure_from_coordinates
from collections import Counter
from geom import Coord, alignment
from dictionaryInspector import DictionaryInspector
from callsys import CallxTB


def build_push_l(coord, at_l, nat1=0):
    # avoiding crossing push atomic pairs
    if len(at_l) == 2:

        at1 = int(at_l[0][0])
        at2 = int(at_l[0][1]) + int(nat1)
        at3 = int(at_l[1][0])
        at4 = int(at_l[1][1]) + int(nat1)
        if Coord(coord).test_if_cross(at1, at2, at3, at4):
            current_atl = [[at1, at4], [at2, at3]]
        else:
            current_atl = [[at1, at2], [at3, at4]]
    else:
        current_atl = [[int(at_l[0][0]), int(at_l[0][1]) + nat1]]

    ## at_l: index starting with 0; push_l: index starting with 1
    push_l = [(str(k[0] + 1),  # at1 (index starts with 1)
               str(k[1] + 1),  # at2 (index starts with 1)
               str(np.round(Coord(coord).get_ats_distance(k[0], k[1]), 2)),  # d_12
               str((Coord(coord).get_rr_distance(k[0]) +
                   Coord(coord).get_rr_distance(k[1])*0.95))) for k in current_atl]  # rr_12

    return push_l


# def standard_scan(scan_wd, coord, at_l, charge, uhf, xtb_setting, mode=2, d=None, tmpat1=None, tmpat2=None,
#                   tmpr=None):
#     os.chdir(scan_wd)
#     level = xtb_setting.get('xtblevel', 'gfn2')
#     alpb = xtb_setting.get('alpb', None)
#     etemp = xtb_setting.get("etemp", str(300))
#
#     write_xyz(coord, "coord")
#
#     if mode == 0:
#         # intraA
#
#         current_at_l = at_l
#         push_l = [(str(k[0] + 1),  # at1 (index starts with 1)
#                    str(k[1] + 1),  # at2 (index starts with 1)
#                    str(np.round(Coord(coord).get_ats_distance(k[0], k[1]), 2)),  # d_12
#                    str(Coord(coord).get_rr_distance(k[0])) +
#                    str(Coord(coord).get_rr_distance(k[1]))) for k in current_at_l]  # rr_12
#
#     elif mode == 1:
#
#         write_xyz(coord, "coord")
#         at1, at2 = at_l[0], at_l[1]
#         # d = searching_setting.get("ddissociation", 2)
#         d = 2
#         d12 = np.round(Coord(coord).get_ats_distance(int(at1), int(at2)), 3)
#
#         push_l = [(str(int(at1) + 1),  # at1
#                    str(int(at2) + 1),  # at2
#                    str(d12),  # d_12
#                    str(d * d12))]  # 2 * d_12
#
#     else:
#         write_xyz(coord, "coord")
#
#     struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
#     scan_setting = {
#         "constrain":
#             {"k": xtb_setting.get("constrain_strength", str(1)),
#              "distance": [str(int(tmpat1) + 1),
#                           str(int(tmpat2) + 1),
#                           str(tmpr * 1.5)] if tmpr is not None else None  # tmpr is a loosed bond
#              },
#         "scan":
#             {"mode": "concerted" if len(at_l) == 2 else None,
#              "distance": push_l if push_l else None,
#              "step": xtb_setting.get("scan_step", str(20))}
#     }
#     rc = struc_callxtb.call_scan(scan_setting)
#     if rc == 0:
#         # scan success
#         return True
#     else:
#         return False
#
# else:
#
#
#
#
# struc_callxtb = CallxTB(level, "coord", charge, uhf, alpb, etemp)
# scan_setting = {
#     "constrain": {
#         "k": xtb_setting.get("constrain_strength", str(1)),
#         "distance": [str(int(tmpat1) + 1),
#                      str(int(tmpat2) + 1),
#                      str(tmpr * 1.5)] if tmpr is not None else None  # tmpr is a loosed bond
#     },
#     "scan": {
#         "mode": "concerted" if len(at_l) == 2 else None,
#         "distance": at_l if at_l else None,  # here at_l is push_l
#         "step": xtb_setting.get("scan_step", str(20))
#     }
# }
# rc = struc_callxtb.call_scan(scan_setting)
# if rc == 0:
#     # scan success
#     return True
# else:
#     return False


def get_endxyz_from_mb(constrain_l, coordinates, ele):
    extra_constraints = {'bond_order_assignment': True,
                         'cycle_detection': True,
                         'set_edges': False,  # False if no bonds should be constrained automatically.
                         'set_angles': True,  # False if no angles should be constrained automatically.
                         'set_dihedrals': False,  # False if no dihedrals should be constrained automatically.
                         'set_repulsion': False,  # False if no coulomb term should be used automatically.
                         'repulsion_charge': 100.0,
                         # Charged used for the Coulomb term in the Force field, every atom-atom pair uses the same charge, only reasonable opt mode (standalone force-field optimization). Defaults to 100.0
                         'constraints': {'bonds': constrain_l}
                         # atoms: list of atoms that define the dihedral, value is the ideal dihedral angle in degrees, atom indexing starts with 1.
                         }
    _ = idealize_structure_from_coordinates(coordinates, ele, input_constraint=extra_constraints)
    coord_l, ele = _[1], _[2]

    coord = pd.DataFrame(coord_l, columns=['X', 'Y', 'Z'])
    coord.insert(0, 'Symbol', ele)
    coord.columns = [''] * len(coord.columns)

    return coord


def evaluate_bo_after_mbff(wd, level, coord, charge, uhf, alpb, etemp, tmp_at_connection, wbolim=0.9):
    # do sp to get wbo
    Trigger = False

    os.chdir(os.path.abspath(wd))
    write_xyz(coord, f"coord")

    rc = CallxTB(level, "coord", charge, uhf, alpb, etemp).callopt()
    if rc==0:
        BO = {}
        wboe = pd.read_csv(os.path.join(wd, "wbo"), header=None, sep='\s+')
        for index in range(wboe.shape[0]):
            dgap = wboe.iloc[index, 2]
            if dgap >= wbolim:
                at_1, at_2 = wboe.iloc[index, 0], wboe.iloc[index, 1]
                value = [at_1, at_2]
                if dgap not in BO.keys():
                    BO[dgap] = [value]
                else:
                    BO[dgap].append(value)
        if tmp_at_connection is not None:
            for item in BO.values():
                for atl1 in item:
                    for atl2 in tmp_at_connection:
                        if Counter(atl1) == Counter(atl2):
                            Trigger = True
        else:
            Trigger = True  # opt sussess

    return Trigger


def get_bo(wbopath,tmp_at_connection):
    wboe = pd.read_csv(os.path.join(wbopath, "wbo"), header=None, sep='\s+')
    BO = {}
    for index in range(wboe.shape[0]):
        dgap = wboe.iloc[index, 2]
        if dgap >= 0.3:
            at_1, at_2 = wboe.iloc[index, 0], wboe.iloc[index, 1]
            value = [at_1, at_2]
            if dgap not in BO.keys():
                BO[dgap] = [value]
            else:
                BO[dgap].append(value)

    bo = 0
    if tmp_at_connection is not None:
        for key,item in BO.items():
            for atl1 in item:
                for atl2 in tmp_at_connection:
                    if Counter(atl1) == Counter(atl2):
                        bo = key
                        continue
    return bo

def get_bl(coord,tmp_at_connection):
    at1 = tmp_at_connection[0][0]-1
    at2 = tmp_at_connection[0][1]-1

    d = np.round(Coord(coord).get_ats_distance(at1, at2), 2)
    return d

