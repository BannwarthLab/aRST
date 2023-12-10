import os
import numpy as np
# from molbar.barcode import idealize_structure_from_file
from molbar.barcode import idealize_structure_from_coordinates
import pandas as pd
from geom import r0


def read_topomol(molfile, nat1=None):
    # V2000
    with open(molfile) as f:
        context = f.readlines()

    nat, nbond = int(context[3].split()[0]), int(context[3].split()[1])
    coordinates = []
    constrain_l = []
    ele = []

    for i in range(nat):
        tmp = context[i + 4].split()[:3]
        coordinates.append([float(k) for k in tmp])
        ele.append(context[i + 4].split()[3])

    for i in range(nbond):
        bond_pair = context[i + nat + 4].split()[:2]
        atom1, atom2 = int(bond_pair[0]), int(bond_pair[1])
        d = r0[ele[atom1 - 1]] + r0[ele[atom2 - 1]]
        tmp_d = {"atoms": {}, "value": {}}
        if nat1 is None:
            tmp_d["atoms"] = [atom1, atom2]
            tmp_d["value"] = d
        else:
            tmp_d["atoms"] = [atom1 + nat1, atom2 + nat1]
            tmp_d["value"] = d
        constrain_l.append(tmp_d)

    return constrain_l, coordinates, ele

def get_endxyz_from_mb(molfile_path,push_l):
    """molfile: xtbtopomol file of aligned coord
       push_l: reactiv atomic list
        """

    constrain_l, coordinates, ele = read_topomol(molfile_path)


    for item in push_l:
        at1 = int(item[0])
        at2 = int(item[1])
        d = float(item[3])

        constrain_l.append({"atoms": [at1, at2], "value":d})


    extra_constraints = {'bond_order_assignment': True,
                         'cycle_detection': True,
                         'set_edges': False,  # False if no bonds should be constrained automatically.
                         'set_angles': False,  # False if no angles should be constrained automatically.
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