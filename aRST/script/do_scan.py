import pandas as pd
from aRST.script.geom import Coord



def build_push_l(coordinate, ele_l, at_l, nat1=0,d=0.98,force_combination=False):
    # avoiding crossing push atomic pairs
    instantce = Coord(coordinate, ele_l)
    if len(at_l) ==2:
        at1 = int(at_l[0][0])
        at2 = int(at_l[0][1]) + int(nat1)
        at3 = int(at_l[1][0])
        at4 = int(at_l[1][1]) + int(nat1)

        if not force_combination and instantce.choose_combination(at1, at2, at3, at4)=="swapped":
            current_atl = [[at1, at4], [at2, at3]]
        else:
            current_atl = [[at1, at2], [at3, at4]]

        # current_atl = [[at1, at2]]
    else:
        current_atl = [[int(at_l[0][0]), int(at_l[0][1]) + nat1]]

    ## at_l: index starting with 0; push_l: index starting with 1
    push_l = [(k[0] + 1,  # at1 (index starts with 1)
               k[1] + 1,  # at2 (index starts with 1)
               instantce.get_ats_distance(k[0], k[1]),  # d_12
               (instantce.rcov_l[k[0]] +instantce.rcov_l[k[1]])*d
               ) for k in current_atl]  # rr_12

    return push_l

def get_endxyz_from_mb(constrain_l, coordinate, ele_l,set_angles=True,custum_constraints=None,repulsion_charge=100.0):
    from molbar.idealize import idealize_structure_from_coordinates
    import numpy as np
    # Define constraints
    if custum_constraints is None:
        extra_constraints = {
            'bond_order_assignment': True,
            'cycle_detection': True,
            'set_edges': False,
            'set_angles': set_angles,
            'set_dihedrals': False,
            'set_repulsion': True,
            'repulsion_charge': repulsion_charge,
            'constraints': {'bonds': constrain_l}
        }
    else:
        extra_constraints = custum_constraints

    # Call the idealization function
    res = idealize_structure_from_coordinates(coordinate, ele_l,
                                              input_constraint=extra_constraints)
    coordinate, ele_l = res[1], res[2]
    return coordinate, ele_l
