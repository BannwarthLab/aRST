import subprocess
import os



import numpy as np


# coordinate = np.array([
#     [-3.0e-14,  0.0e+00, 2.04974830768590],
#     [ 5.0e-14, -3.0e-14, 2.98421207310016],
#     [-2.0e-14,  2.0e-14, 3.91900080449870],
# ], dtype=np.float64)
#
# ele_l = np.array(['H', 'O', 'H'])
#
#
# topo_l = [
#     {
#         "atoms": [1, 2],
#         "value": 0.93446376541426
#     },
#     {
#         "atoms": [2, 3],
#         "value": 0.93478873139854
#     }
# ]
# from molbar.idealize import idealize_structure_from_coordinates
# # Define constraints
# extra_constraints = {
#     'bond_order_assignment': True,
#     'cycle_detection': True,
#     'set_edges': False,
#     'set_angles': True,
#     'set_dihedrals': False,
#     'set_repulsion': True,
#     'repulsion_charge': 100,
#     'constraints': {'bonds': topo_l}
# }
# print(extra_constraints)
# # Call the idealization function
# res = idealize_structure_from_coordinates(coordinate,ele_l,
#                                           input_constraint=extra_constraints)
# coordinate, ele_l = res[1], res[2]
# print(coordinate, ele_l)



