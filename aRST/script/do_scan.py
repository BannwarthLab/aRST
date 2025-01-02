import importlib
import os
import sys
import numpy as np
import pandas as pd
import importlib
from aRST.script.toolbox import  write_xyz
from molbar.barcode import idealize_structure_from_coordinates
import copy
from collections import Counter
from aRST.script.geom import Coord
from aRST.script.callsys import CallxTB


def build_push_l(coord, at_l, nat1=0,d=0.98):
    # avoiding crossing push atomic pairs
    if len(at_l) ==2:

        at1 = int(at_l[0][0])
        at2 = int(at_l[0][1]) + int(nat1)
        at3 = int(at_l[1][0])
        at4 = int(at_l[1][1]) + int(nat1)
        if Coord(coord).test_if_cross(at1, at2, at3, at4):
            current_atl = [[at1, at4], [at2, at3]]
        else:
            current_atl = [[at1, at2], [at3, at4]]
        # current_atl = [[at1, at2]]
    else:
        current_atl = [[int(at_l[0][0]), int(at_l[0][1]) + nat1]]

    ## at_l: index starting with 0; push_l: index starting with 1
    push_l = [(k[0] + 1,  # at1 (index starts with 1)
               k[1] + 1,  # at2 (index starts with 1)
               Coord(coord).get_ats_distance(k[0], k[1]),  # d_12
               (Coord(coord).get_rr_distance(k[0]) +
                   Coord(coord).get_rr_distance(k[1])*d)) for k in current_atl]  # rr_12

    return push_l

def get_endxyz_from_mb(constrain_l, coordinates, ele,set_angles=True):
    # Define constraints
    extra_constraints = {
        'bond_order_assignment': True,
        'cycle_detection': True,
        'set_edges': False,
        'set_angles': set_angles,
        'set_dihedrals': False,
        'set_repulsion': False,
        'repulsion_charge': 100.0,
        'constraints': {'bonds': constrain_l}
    }

    # Call the idealization function
    res = idealize_structure_from_coordinates(coordinates, ele, input_constraint=extra_constraints)
    coord_l, ele = res[1], res[2]

    # Format results into a dataframe
    coord = pd.DataFrame(coord_l, columns=['X', 'Y', 'Z'])
    coord.insert(0, 'Symbol', ele)
    coord.columns = [''] * len(coord.columns)
    return coord


# def evaluate_bo_after_mbff(wd, level, coord, charge, uhf, alpb, etemp, tmp_at_connection, wbolim=0.9):
#     # do sp to get wbo
#     Trigger = False
#
#     os.chdir(os.path.abspath(wd))
#     write_xyz(coord, f"coord")
#
#     rc = CallxTB(level, "coord", charge, uhf, alpb, etemp).callopt()
#     if rc==0:
#         BO = {}
#         wboe = pd.read_csv(os.path.join(wd, "wbo"), header=None, sep='\s+')
#         for index in range(wboe.shape[0]):
#             dgap = wboe.iloc[index, 2]
#             if dgap >= wbolim:
#                 at_1, at_2 = wboe.iloc[index, 0], wboe.iloc[index, 1]
#                 value = [at_1, at_2]
#                 if dgap not in BO.keys():
#                     BO[dgap] = [value]
#                 else:
#                     BO[dgap].append(value)
#         if tmp_at_connection is not None:
#             for item in BO.values():
#                 for atl1 in item:
#                     for atl2 in tmp_at_connection:
#                         if Counter(atl1) == Counter(atl2):
#                             Trigger = True
#         else:
#             Trigger = True  # opt sussess
#
#     return Trigger


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

# def get_bl(coord,tmp_at_connection):
#     at1 = tmp_at_connection[0][0]-1
#     at2 = tmp_at_connection[0][1]-1
#
#     d = np.round(Coord(coord).get_ats_distance(at1, at2), 2)
#     return d

