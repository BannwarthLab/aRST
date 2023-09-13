import sys
import pandas as pd
from spr import *
from geom import alignment
from geom import Coord
import random
from operation import creat_folder, write_xyz, boltzmanp
from callsys import call_xtb_scan


def push(reactants, rec, wd, step=0, **kwargs):
    if kwargs.get('kwargs'):
        kwargs = kwargs.get('kwargs')
    alpb = kwargs.get('alpb', None)
    reactor = kwargs.get('reactor_like', False)
    mode = kwargs.get('platform', 0)
    constrain_strength = kwargs.get('constrain_strength', str(1))
    reactants_dict = {}

    # radical_trigger = False
    # for name in reactants:
    #     if rec[name]['m'] != 1:
    #         radical_trigger = True

    # get all association possibilities
    push_dict_1, push_dict_2 = spr_inter_new(reactants=reactants, rec=rec)

    if reactor:
        # random sample one association
        gap_l = cut_push_dict(push_dict_1, push_dict_2, reactants, 4)

        if gap_l:
            # add boltzman sampling to evaluate min_gap
            boltzp_l = [boltzmanp(k) for k in gap_l]
            boltzp_l = [k / sum(boltzp_l) for k in boltzp_l]
            z = random.random()
            if boltzp_l[0] > z:
                gap = sorted(list(push_dict_1.keys()))[0]
            else:
                gap = None
            if gap is not None:
                # detemine active molecules and atoms
                struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, gap, reactants)

                # get molecular info
                coord_1 = rec[struc_1]['coord']
                nat_1 = coord_1.shape[0]
                coord_2 = rec[struc_2]['coord']

                # write reactants xyz Not used, only for checking
                reactants_dict[f'{str(step + 1)}_1'] = [struc_1, struc_2]
                creat_folder(f'{wd}/{str(step + 1)}_1/')
                write_xyz(coord_1, "coord1")
                write_xyz(coord_2, "coord2")

                # write scan initial xyz
                coord_merge, rr = alignment(coord_1, coord_2, at_1, at_2)
                creat_folder(f'{wd}/{str(step + 1)}_1/scan')
                write_xyz(coord_merge, 'coord')
                with open("gap",'w') as f:
                    f.write(str(gap))

                # submit scan job
                charg = int(rec[struc_1]['charge'] + rec[struc_2]['charge'])
                if int(rec[struc_1]['m']) == 1 and int(rec[struc_2]['m']) == 1:
                    uhf = 0
                else:
                    uhf = int(rec[struc_1]['m']) + int(rec[struc_2]['m']) - 2

                if len(at_l) > 1:
                    at_3, at_4 = int(at_l[1][0]), int(at_l[1][1])
                    r3 = Coord(eval('coord_merge')).get_punkt_coord(at_3)
                    r4 = Coord(eval('coord_merge')).get_punkt_coord(at_4 + nat_1)
                    r = np.linalg.norm(r3 - r4)
                    call_xtb_scan(f'{wd}/{str(step + 1)}_1/scan', charg, uhf, alpb=alpb, mode=mode,
                                  kw={'at3': at_3, 'at4': at_4 + nat_1, 'concerted_trigger': True,
                                      'at1': at_1, 'at2': at_2 + nat_1, 'd1': rr * 5, 'd2': rr,
                                      'd_34_start': r, 'd_34_end': rr, 'constrain_strength': constrain_strength})
                else:
                    call_xtb_scan(f'{wd}/{str(step + 1)}_1/scan', charg, uhf, alpb=alpb, mode=mode,
                                  kw={'concerted_trigger': False, 'at1': at_1, 'at2': at_2 + nat_1,
                                      'd1': rr * 5, 'd2': rr, 'constrain_strength': constrain_strength})
        else:
            print("gap_l wrong: ", reactants, push_dict_1)
    else:
        # full search
        candidates = kwargs.get('intercands', 1)
        gap_l = cut_push_dict(push_dict_1, push_dict_2, reactants, candidates)
        count = 1
        if gap_l:
            for i in gap_l:
                # get reactants name and active atoms
                struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, i, reactants)

                # get correspongding molecular info
                coord_1 = rec[struc_1]['coord']
                nat_1 = coord_1.shape[0]
                coord_2 = rec[struc_2]['coord']
                nreact = struc_1.split('/')[0] + str(count)
                reactants_dict[nreact] = [struc_1, struc_2]

                print("look hier!! ", nreact, struc_1, struc_2, at_1, at_2, i, coord_1.iloc[at_1 - 1, 0],
                      coord_2.iloc[at_2 - 1, 0])

                # do alignment and write scan starting .xyz
                coord_merge, rr = alignment(coord_1, coord_2, at_1, at_2)
                creat_folder(f'{wd}/{nreact}/scan')
                write_xyz(coord_merge, 'coord')

                with open("gap",'w') as f:
                    f.write(str(i))

                # write candidates posiibility for checking up
                # df = pd.DataFrame.from_dict(push_dict_1, orient='index')
                # df.to_csv(f'{wd}/{nreact}/scan/candidates.csv', header=False)

                charg = int(rec[struc_1]['charge'] + rec[struc_2]['charge'])
                if int(rec[struc_1]['m']) == 1 and int(rec[struc_2]['m']) == 1:
                    uhf = 0
                else:
                    uhf = int(rec[struc_1]['m']) + int(rec[struc_2]['m']) - 2
                    # tmp setting
                    if uhf == 2:
                        uhf = 0

                if len(at_l) > 1:
                    at_3, at_4 = int(at_l[1][0]), int(at_l[1][1])
                    r3 = Coord(eval(f'coord_merge')).get_punkt_coord(at_3)
                    r4 = Coord(eval(f'coord_merge')).get_punkt_coord(at_4 + nat_1)
                    r = np.linalg.norm(r3 - r4)

                    call_xtb_scan(f'{wd}/{nreact}/scan', charg, uhf, alpb=alpb, mode=mode,
                                  kw={'at3': at_3, 'at4': at_4 + nat_1, 'concerted_trigger': True,
                                      'at1': at_1, 'at2': at_2 + nat_1, 'd1': rr * 5, 'd2': rr,
                                      'd_34_start': r, 'd_34_end': rr, 'constrain_strength': constrain_strength,
                                      'etmp':kwargs.get('etmp', str(300))})
                else:
                    call_xtb_scan(f'{wd}/{nreact}/scan', charg, uhf, alpb=alpb, mode=mode,
                                  kw={'concerted_trigger': False,
                                      'at1': at_1, 'at2': at_2 + nat_1, 'd1': rr * 5, 'd2': rr,
                                      'constrain_strength': constrain_strength,'etmp':kwargs.get('etmp', str(300))})

                count += 1
    return reactants_dict


def pull(reactants, rec, wd, step=0, **kwargs):
    if kwargs.get('kwargs'):
        kwargs = kwargs.get('kwargs')
    alpb = kwargs.get('alpb', None)
    reactor = kwargs.get('reactor_like', False)
    hydrogen = kwargs.get('hydrogen', True)
    mode = kwargs.get('platform', 0)
    constrain_strength = kwargs.get('constrain_strength', str(1))
    discands = kwargs.get('discands', 1)

    reactants_dict = {}
    if not hydrogen:
        pull_dict = spr_intra_new(reactants, rec, hydrogen=False)
    else:
        pull_dict = spr_intra_new(reactants, rec)
    if pull_dict:
        # # add boltzmann sampling
        # gap_l = sorted(list(pull_dict.keys()))[:5]
        # boltzp_l = [boltzmanp(k) for k in gap_l]
        # boltzp_l = [k / sum(boltzp_l) for k in boltzp_l]
        # z = random.random()
        # if boltzp_l[0] > z:
        #     gap = min(list(pull_dict.keys()))
        # else:
        #     gap = None
        # gap = min(list(pull_dict.keys()))
        gap_l = list(sorted(list(pull_dict.keys())))[:discands]
        if gap_l:
            for i in range(len(gap_l)):
                gap = gap_l[i]
                # detemine active molecule and atoms
                # min_gap = min(list(pull_dict.keys()))
                dis_struc = pull_dict[gap][0]
                dis_at_1 = pull_dict[gap][1]
                dis_at_2 = pull_dict[gap][2]
                coord = rec[dis_struc]['coord']
                C = Coord(coord)
                r1 = C.get_punkt_coord(dis_at_1)
                r2 = C.get_punkt_coord(dis_at_2)
                r = np.linalg.norm(r1 - r2)
                charg = rec[dis_struc]['charge']
                m = rec[dis_struc]['m']
                uhf = m - 1
                if reactor and i == 0:
                    # get reactant name and active atoms
                    reactants_dict[f'{str(step + 1)}_2'] = [dis_struc]

                    # write reactant xyz
                    creat_folder(f'{wd}/{str(step + 1)}_2')
                    dis_struc_tmp = f'{dis_struc}'.replace('/', '_')
                    write_xyz(coord, f'{dis_struc_tmp}')  # Not used, only for checking

                    # submit scan job
                    creat_folder(f'{wd}/{str(step + 1)}_2/scan')
                    write_xyz(coord, 'coord')

                    call_xtb_scan(f'{wd}/{str(step + 1)}_2/scan', charg, uhf=uhf, alpb=alpb, mode=mode,
                                  kw={'concerted_trigger': False, 'at1': dis_at_1, 'at2': dis_at_2, 'd1':
                                      r, 'd2': r * 2, 'constrain_strength': constrain_strength,
                                      'etmp':kwargs.get('etmp', str(300))})

                else:
                    nreact_count = str(kwargs.get('association_count', 5) + 1+i)
                    nreact = dis_struc.split('/')[0] + nreact_count
                    reactants_dict[nreact] = [dis_struc]
                    creat_folder(f'{wd}/{nreact}/scan')
                    write_xyz(coord, 'coord')
                    call_xtb_scan(f'{wd}/{nreact}/scan', charg, uhf=uhf, alpb=alpb, mode=mode,
                                  kw={'concerted_trigger': False,
                                      'at1': dis_at_1, 'at2': dis_at_2, 'd1': r, 'd2': r * 1.5,
                                      'constrain_strength': constrain_strength,'etmp':kwargs.get('etmp', str(300))})
    return reactants_dict


def assign_asite(react_dict1, react_dict2, gap_i, reactants):
    # detemine active molecules and atoms and pass on to the next function (only for association)
    # struc: reactants
    # at: active atoms
    # at_l: possibile equivalent atomic pairs for reaction sites
    if gap_i in react_dict1.keys():
        react_dict = react_dict1
    else:
        react_dict = react_dict2
    struc_1_nr, struc_2_nr = ast.literal_eval(react_dict[gap_i][0])[0], ast.literal_eval(react_dict[gap_i][0])[1]
    struc_1, struc_2 = reactants[int(struc_1_nr) - 1], reactants[int(struc_2_nr) - 1]
    at_l = ast.literal_eval(react_dict[gap_i][1])
    at_1, at_2 = at_l[0][0], at_l[0][1]
    at_1, at_2 = int(at_1), int(at_2)

    return struc_1, struc_2, at_1, at_2, at_l


def cut_push_dict(push_dict_1, push_dict_2, reactants, candidates):
    if candidates > 1:
        gap_l1 = sorted(list(push_dict_1.keys()))
        gap_l1 = [gap_l1[0]] + list(sorted([gap_l1[i] for i in range(len(gap_l1))
                                            if i > 0 and gap_l1[i] - gap_l1[i - 1] > 1e-3]))[:int(candidates / 2) - 1]
        tmp_l = []
        for k in gap_l1:
            struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, k, reactants)
            tmp_l.append([struc_1, struc_2, at_1, at_2])


        gap_l2 = sorted(list(push_dict_2.keys()))
        tmp_l3 = []
        struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, gap_l2[0], reactants)
        if [struc_1, struc_2, at_1, at_2] not in tmp_l and [struc_2, struc_1, at_2, at_1] not in tmp_l:
            tmp_l3.append(gap_l2[0])
        for k in range(1, len(gap_l2) - 1):
            if gap_l2[k] - gap_l2[k - 1] > 1e-3:
                struc_1, struc_2, at_1, at_2, at_l = assign_asite(push_dict_1, push_dict_2, gap_l2[k], reactants)
                if [struc_1, struc_2, at_1, at_2] not in tmp_l and [struc_2, struc_1, at_2, at_1] not in tmp_l:
                    tmp_l3.append(gap_l2[k])
        tmp_l3 = tmp_l3[:int(candidates / 2)]
        gap_l = gap_l1 + tmp_l3
    else:
        gap_l = [[sorted(list(push_dict_1.keys()))][0]]
    return gap_l
