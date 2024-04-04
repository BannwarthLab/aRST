# # from molbar.barcode import idealize_structure_from_file
# from molbar.barcode import idealize_structure_from_coordinates
# import pandas as pd
# from toolbox import read_topomol
#
#
# def get_endxyz_from_mb_old(molfile_path, push_l, tmpr_l):
#     """molfile: xtbtopomol file of aligned coord
#        push_l: reactiv atomic list
#         """
#
#     constrain_l, coordinates, ele = read_topomol(molfile_path)
#
#     for item in push_l:
#         at1 = int(item[0])
#         at2 = int(item[1])
#         d = float(item[3])
#
#         constrain_l.append({"atoms": [at1, at2], "value": d})
#         # constrain_bond_l.append({"atoms": [at1, at2], "value":d})
#
#     if tmpr_l is not None:  # loose bond
#         for item in tmpr_l:
#             at1 = int(item[0])
#             at2 = int(item[1])
#             d = float(item[2])
#             constrain_l.append({"atoms": [at1, at2], "value": d})
#
#     with open("constrain_l", "w") as f:
#         for i in constrain_l:
#             f.write(str(i) + "\n")
#
#     extra_constraints = {'bond_order_assignment': True,
#                          'cycle_detection': True,
#                          'set_edges': True,  # False if no bonds should be constrained automatically.
#                          'set_angles': True,  # False if no angles should be constrained automatically.
#                          'set_dihedrals': True,  # False if no dihedrals should be constrained automatically.
#                          'set_repulsion': True,  # False if no coulomb term should be used automatically.
#                          'repulsion_charge': 100.0,
#                          # Charged used for the Coulomb term in the Force field, every atom-atom pair uses the same charge, only reasonable opt mode (standalone force-field optimization). Defaults to 100.0
#                          'constraints': {'bonds': constrain_l}
#                          # atoms: list of atoms that define the dihedral, value is the ideal dihedral angle in degrees, atom indexing starts with 1.
#                          }
#     _ = idealize_structure_from_coordinates(coordinates, ele, input_constraint=extra_constraints)
#     coord_l, ele = _[1], _[2]
#
#     coord = pd.DataFrame(coord_l, columns=['X', 'Y', 'Z'])
#     coord.insert(0, 'Symbol', ele)
#     coord.columns = [''] * len(coord.columns)
#
#     return coord
#
#
#
#
#
#
#
#
#
#
#     coord = struc_inspector.get_value(f"{strucid}.struc_info.general.coord")
#
#     write_xyz(coord, "coord")
#     charge = struc_inspector.get_value(f"{strucid}.struc_info.general.charge")
#     uhf = struc_inspector.get_value(f"{strucid}.struc_info.general.multiplicity") - 1
#
#     at1, at2, at3, at4 = int(at_l[0][0]), int(at_l[0][1]), \
#                          int(at_l[1][0]), int(at_l[1][1])
#     if len(at_l) == 2 and Coord(coord).test_if_cross(at1, at2, at3, at4):
#         current_at_l = [[at1, at4], [at2, at3]]
#     elif len(at_l) == 2:
#         current_at_l = [[at1, at2], [at3, at4]]
#     else:
#         current_at_l = [[at1, at2]]
#
#
# push_l = [(str(int(k[0]) + 1),  # at1 (index starts with 1)
#            str(int(k[1]) + 1),  # at2 (index starts with 1)
#            str(np.round(Coord(coord).get_ats_distance(int(k[0]), int(k[1])), 3)),  # d_12
#            str(Coord(coord).get_rr_distance(int(k[0])) +
#                Coord(coord).get_rr_distance(int(k[1])))  # rr_12
#            ) for k in current_at_l]
#
# tmp_at_connection = []
# for item in push_l:
#     at1 = int(item[0])
#     at2 = int(item[1])
#     d = float(item[3])
#
#     constrain_l.append({"atoms": [at1, at2], "value": d})
#
#     tmp_at_connection.append([at1, at2])
#
# if tmp_l is not None:  # loose bond
#     for item in tmp_l:
#         at1 = int(item[0])
#         at2 = int(item[1])
#         d = float(item[2])
#         constrain_l.append({"atoms": [at1, at2], "value": d})
#
# with open("constrain_l", "w") as f:
#     for i in constrain_l:
#         f.write(str(i) + "\n")
# extra_constraints = {'bond_order_assignment': True,
#                      'cycle_detection': True,
#                      'set_edges': False,  # False if no bonds should be constrained automatically.
#                      'set_angles': True,  # False if no angles should be constrained automatically.
#                      'set_dihedrals': False,  # False if no dihedrals should be constrained automatically.
#                      'set_repulsion': False,  # False if no coulomb term should be used automatically.
#                      'repulsion_charge': 100.0,
#                      # Charged used for the Coulomb term in the Force field, every atom-atom pair uses the same charge, only reasonable opt mode (standalone force-field optimization). Defaults to 100.0
#                      'constraints': {'bonds': constrain_l}
#                      # atoms: list of atoms that define the dihedral, value is the ideal dihedral angle in degrees, atom indexing starts with 1.
#                      }
# coordinates = []
# for line in range(coord.shape[0]):
#     coordinates.append(coord.loc[line, 1:4].tolist())
#
# _ = idealize_structure_from_coordinates(coordinates, ele, input_constraint=extra_constraints)
# coord_l, ele = _[1], _[2]
#
# endxyz = pd.DataFrame(coord_l, columns=['X', 'Y', 'Z'])
# endxyz.insert(0, 'Symbol', ele)
# endxyz.columns = [''] * len(endxyz.columns)
# write_xyz(endxyz, "endxyz")
