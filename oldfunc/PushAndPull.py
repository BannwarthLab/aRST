# import copy
# import sys
# import pandas as pd
# from spr import *
# from geom import alignment
# from geom import Coord
# import random
# from operation import create_folder, write_xyz, boltzmanp
#
# from collections import Counter
#
#
# def association_pairs(reactants, rec, **kwargs):
#     reactor = kwargs.get('reactor_like', False)
#
#     # get all association possibilities
#     push_dict_1, push_dict_2, hoaoa_df, hoaob_df, luaoa_df, luaob_df = spr_inter_new(reactants=reactants, rec=rec)
#     output_d = {}
#     if reactor:
#         # random sample one association
#         gap_l = cut_push_dict(push_dict_1, push_dict_2, reactants, 4)
#
#         if gap_l:
#             # add boltzman sampling to evaluate min_gap
#             boltzp_l = [boltzmanp(k) for k in gap_l]
#             boltzp_l = [k / sum(boltzp_l) for k in boltzp_l]
#             z = random.random()
#             if boltzp_l[0] > z:
#                 gap = sorted(list(push_dict_1.keys()))[0]
#             else:
#                 gap = None
#             if gap is not None:
#                 # detemine active molecules and atoms
#                 struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, gap, reactants)
#                 output_d[gap] = [struc_1, struc_2, at_1, at_2, at_l]
#     else:
#         # full search
#         candidates = kwargs.get('intermolecular_pair', 1)
#
#         gap_l = cut_push_dict(push_dict_1, push_dict_2, reactants, len(push_dict_1))
#
#         # tmp seting for all reaction possibilities (25.10.2023)
#         if candidates == 'all':
#             for i in gap_l:
#                 struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, i, reactants)
#                 output_d[i] = struc_1, struc_2, at_1, at_2, at_l
#         elif candidates != 0:
#             count_2 = 0  # each 2 reactants in the system
#             tmp_l = []
#             if gap_l:
#                 loop_count = 0
#                 while loop_count < len(gap_l) and count_2 < Cn_2(len(reactants)):
#                     i = gap_l[loop_count]
#                     loop_count += 1
#                     # get reactants name and active atoms
#                     struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, i, reactants)
#                     # count how many candidates have been selected for struc1 and struct2
#                     tmp_c = 0
#                     for j in tmp_l:
#                         if Counter([struc_1, struc_2]) == Counter(j):
#                             tmp_c += 1
#                     tmp_l.append([struc_1, struc_2])
#                     if tmp_c < candidates:
#                         output_d[i] = struc_1, struc_2, at_1, at_2, at_l
#                     elif tmp_c == candidates:
#                         count_2 += 1
#     return output_d, hoaoa_df, hoaob_df, luaoa_df, luaob_df
#
#
# def push(association_d, rec, wd, reactor, step=0, **xtb_setting):
#     reactants_dict = {}
#     count = 1
#     for gap in list(association_d.keys()):
#         struc_1, struc_2, at_1, at_2, at_l = association_d[gap]
#
#         # get molecular info
#         coord_1 = rec[struc_1]['coord']
#         nat_1 = coord_1.shape[0]
#         coord_2 = rec[struc_2]['coord']
#
#         if reactor:  # todo _ab?
#             create_folder(f'{wd}/{str(step + 1)}_1/scan')
#             filename = f'{wd}/{str(step + 1)}_1/scan'
#             reactants_dict[f'{str(step + 1)}_1'] = [struc_1, struc_2]  # write reactants xyz Not used, only for checking
#
#         else:
#             nreact = struc_1.split('/')[0] + '_' + str(count)
#             create_folder(f'{wd}/{nreact}/scan')
#             filename = f'{wd}/{nreact}/scan'
#             reactants_dict[nreact] = [struc_1, struc_2]  # write reactants xyz Not used, only for checking
#         # write scan initial xyz
#         # test:
#         write_xyz(coord_1, "coord1")
#         write_xyz(coord_2, "coord2")
#
#         coord_merge, rr = alignment(coord_1, coord_2, at_1, at_2)
#         write_xyz(coord_merge, 'coord')
#         with open("gap", 'w') as f:
#             f.write(str(gap))
#
#         # submit scan job
#         charg = int(rec[struc_1]['charge'] + rec[struc_2]['charge'])
#         if int(rec[struc_1]['m']) == 1 and int(rec[struc_2]['m']) == 1:
#             uhf = 0
#         elif int(rec[struc_1]['m']) == 2 and int(rec[struc_2]['m']) == 2:
#             uhf = 0
#         else:
#             uhf = int(rec[struc_1]['m']) + int(rec[struc_2]['m']) - 2
#
#         if len(at_l) > 1:
#             at_3, at_4 = int(at_l[1][0]), int(at_l[1][1])
#             r3 = Coord(eval('coord_merge'), at=None).get_punkt_coord(at_3)
#             r4 = Coord(eval('coord_merge'), at=None).get_punkt_coord(at_4 + nat_1)
#             r = np.round(np.linalg.norm(r3 - r4), 3)
#             rr = np.round(rr, 3)
#             tmp_d = copy.deepcopy(xtb_setting)
#             tmp_d.update({'at3': at_3, 'at4': at_4 + nat_1, 'concerted_trigger': True,
#                           'at1': at_1, 'at2': at_2 + nat_1, 'd1': rr * 5, 'd2': rr * 0.9,
#                           'd_34_start': r, 'd_34_end': rr * 0.9})
#             call_xtb_scan(filename, charg, uhf, **tmp_d)
#         else:
#             tmp_d = copy.deepcopy(xtb_setting)
#             tmp_d.update({'concerted_trigger': False, 'at1': at_1, 'at2': at_2 + nat_1,
#                           'd1': rr * 5, 'd2': rr * 0.9})
#             call_xtb_scan(filename, charg, uhf, **tmp_d)
#         count += 1
#     return reactants_dict
#
#
# def dissociation_pair(reactants, rec, **kwargs):
#     hydrogen = kwargs.get('hydrogen', True)
#     reactor = kwargs.get('reactor_like', False)
#     candidates = kwargs.get('dissociation', 1)
#     dissociation_d = {}
#     mode = kwargs.get('dissociation_mode', 0)
#     if not hydrogen:
#         pull_dict = spr_intra_new(reactants, rec, mode=mode, hydrogen=False)
#     else:
#         pull_dict = spr_intra_new(reactants, rec)
#     if reactor:
#         gap = sorted(list(pull_dict.keys()))[0]
#         dis_struc = pull_dict[gap][0]
#         dis_at_1 = pull_dict[gap][1]
#         dis_at_2 = pull_dict[gap][2]
#         dissociation_d[gap] = [dis_struc, dis_at_1, dis_at_2]
#     elif candidates != 0:
#         if candidates == 'all':  # tmp seting for all reaction possibilities
#             candidates = len(pull_dict.keys())
#         sortedgap = sorted(list(pull_dict.keys()))
#         gap_l = []
#         prev_key = None
#         for key in sortedgap:
#             if prev_key is None or key - prev_key > 1e-3:
#                 gap_l.append(key)
#             if len(gap_l) == candidates:
#                 break
#         # gap_l = list(sorted(list(pull_dict.keys())))[:candidates]
#         if gap_l:
#             for i in range(len(gap_l)):
#                 gap = gap_l[i]
#                 # detemine active molecule and atoms
#                 # min_gap = min(list(pull_dict.keys()))
#                 dis_struc = pull_dict[gap][0]
#                 dis_at_1 = pull_dict[gap][1]
#                 dis_at_2 = pull_dict[gap][2]
#                 dissociation_d[gap] = dis_struc, dis_at_1, dis_at_2
#     return dissociation_d
#
#
# def pull(dissociation_d, rec, wd, reactor, ninter, step=0, **xtb_setting):
#     reactants_dict = {}
#     count = 1
#     for gap in list(dissociation_d.keys()):
#         dis_struc, dis_at_1, dis_at_2 = dissociation_d[gap]
#         coord = rec[dis_struc]['coord']
#         C = Coord(coord, at=None)
#         r1 = C.get_punkt_coord(dis_at_1)
#         r2 = C.get_punkt_coord(dis_at_2)
#         r = np.round(np.linalg.norm(r1 - r2),3)
#         charg = rec[dis_struc]['charge']
#         m = rec[dis_struc]['m']
#         uhf = m - 1
#         if reactor:
#             create_folder(f'{wd}/{str(step + 1)}_2/scan')
#             filename = f'{wd}/{str(step + 1)}_2/scan'
#             reactants_dict[f'{str(step + 1)}_2'] = [dis_struc]
#
#         else:
#             nreact_count = str(count + ninter)
#             nreact = dis_struc.split('/')[0] + "_" + nreact_count
#             reactants_dict[nreact] = [dis_struc]
#             create_folder(f'{wd}/{nreact}/scan')
#             filename = f'{wd}/{nreact}/scan'
#
#         dis_struc_tmp = f'{dis_struc}'.replace('/', '.')
#         write_xyz(coord, f'{dis_struc_tmp}')
#         write_xyz(coord, 'coord')
#         d = xtb_setting.get("dissociation_distance", 2)
#
#         tmp_d = copy.deepcopy(xtb_setting)
#         tmp_d.update({'concerted_trigger': False,
#                       'at1': dis_at_1, 'at2': dis_at_2, 'd1': r, 'd2': r * d})
#         call_xtb_scan(filename, charg, uhf, **tmp_d)
#         count += 1
#
#     return reactants_dict
#
#
# def assign_asite(react_dict1, react_dict2, gap_i, reactants):
#     # detemine active molecules and atoms and pass on to the next function (only for association)
#     # struc: reactants
#     # at: active atoms
#     # at_l: possibile equivalent atomic pairs for reaction sites
#     if gap_i in react_dict1.keys():
#         react_dict = react_dict1
#     else:
#         react_dict = react_dict2
#     struc_1_nr, struc_2_nr = ast.literal_eval(react_dict[gap_i][0])[0], ast.literal_eval(react_dict[gap_i][0])[1]
#     struc_1, struc_2 = reactants[int(struc_1_nr) - 1], reactants[int(struc_2_nr) - 1]
#     at_l = ast.literal_eval(react_dict[gap_i][1])
#     at_1, at_2 = at_l[0][0], at_l[0][1]
#     at_1, at_2 = int(at_1), int(at_2)
#
#     return struc_1, struc_2, at_1, at_2, at_l
#
#
# def cut_push_dict(push_dict_1, push_dict_2, reactants, candidates=None):
#     if candidates > 1:
#         gap_l1 = sorted(list(push_dict_1.keys()))
#         gap_l1 = [gap_l1[0]] + list(sorted([gap_l1[i] for i in range(len(gap_l1))
#                                             if i > 0 and gap_l1[i] - gap_l1[i - 1] > 1e-3]))[:int(candidates / 2) - 1]
#
#         tmp_l = []
#         for k in gap_l1:
#             struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, k, reactants)
#             tmp_l.append([struc_1, struc_2, at_1, at_2])
#
#         gap_l2 = sorted(list(push_dict_2.keys()))
#
#         tmp_l3 = []
#         struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, gap_l2[0], reactants)
#         if [struc_1, struc_2, at_1, at_2] not in tmp_l and [struc_2, struc_1, at_2, at_1] not in tmp_l:
#             tmp_l3.append(gap_l2[0])
#         for k in range(1, len(gap_l2)):
#             if gap_l2[k] - gap_l2[k - 1] > 1e-3:
#                 struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, gap_l2[k], reactants)
#                 if [struc_1, struc_2, at_1, at_2] not in tmp_l and [struc_2, struc_1, at_2, at_1] not in tmp_l:
#                     tmp_l3.append(gap_l2[k])
#
#         tmp_l3 = tmp_l3[:int(candidates / 2)]
#         gap_l = gap_l1 + tmp_l3
#     elif push_dict_1.keys():
#         gap_l = [sorted(list(push_dict_1.keys()))[0]]
#     else:
#         gap_l = []
#     return gap_l
#
#
# def Cn_2(n):
#     out = 1
#     for i in range(1, n + 1):
#         out *= i
#     out /= 2
#     return int(out)
