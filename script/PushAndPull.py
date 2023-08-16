from spr import *
from geom import alignment
from geom import Coord
import random
# from callsys import creat_folder, write_xyz, submit_scan
from submit_tmp import creat_folder, write_xyz, submit_scan


def push(reactants, rec, wd, step=0, **kwargs):
    kwargs = kwargs.get('kwargs')
    alpb = kwargs.get('alpb', None)
    reactor = kwargs.get('reactor_like', False)
    reactants_dict = {}

    radical_trigger = False
    for name in reactants:
        if rec[name]['m'] != 1:
            radical_trigger = True

    # get all association possibilities
    push_dict_1, push_dict_2 = spr_inter_new(reactants=reactants, wd_pre=f'{wd}')

    if reactor:
        # random sample one association
        if not radical_trigger:
            gap_l = sorted(list(push_dict_1.keys()))[:5]
        else:
            gap_l = sorted(list(push_dict_1.keys()))[:2] + sorted(list(push_dict_2.keys()))[:2]
            push_dict_1.update(push_dict_2)

        if gap_l:
            #gap = random.sample(gap_l, 1)[0]  # random pick
            gap = sorted(list(push_dict_1.keys()))[0]


            # detemine active molecules and atoms
            struc_1_nr, struc_2_nr = ast.literal_eval(push_dict_1[gap][0])[0], ast.literal_eval(push_dict_1[gap][0])[1]
            struc_1, struc_2 = reactants[int(struc_1_nr) - 1], reactants[int(struc_2_nr) - 1]
            at_l = ast.literal_eval(push_dict_1[gap][1])
            at_1, at_2 = at_l[0][0], at_l[0][1]
            at_1, at_2 = int(at_1), int(at_2)
            coord_1 = rec[struc_1]['coord']
            nat_1 = coord_1.shape[0]
            coord_2 = rec[struc_2]['coord']

            # write reactants xyz
            reactants_dict[f'{str(step + 1)}_1'] = struc_1, struc_2
            # creat_folder(f'{wd}/{str(step + 1)}_1/')
            # for i in range(1, 3):
            #     coord_name = f'coord_{i}'
            #     struc_name = f'struc_{i}'
            #     write_xyz(eval('coord_name'), eval('struc_name'))

            # write scan initial xyz
            coord_merge, rr = alignment(coord_1, coord_2, at_1, at_2)
            creat_folder(f'{wd}/{str(step + 1)}_1/scan')
            write_xyz(coord_merge, 'coord')

            # submit scan job
            charg = int(rec[struc_1]['charge'] + rec[struc_2]['charge'])
            if len(at_l) > 1:
                at_3, at_4 = int(at_l[1][0]), int(at_l[1][1])
                r3 = Coord(eval('coord_merge')).get_punkt_coord(at_3)
                r4 = Coord(eval('coord_merge')).get_punkt_coord(at_4 + nat_1)
                r = np.linalg.norm(r3 - r4)

                submit_scan(f'{wd}/{str(step + 1)}_1/scan', charg, alpb,
                            kw={'at3': at_3, 'at4': at_4 + nat_1, 'concerted_trigger': True,
                                'at1': at_1, 'at2': at_2 + nat_1, 'd1': rr * 5, 'd2': rr,
                                'd_34_start': r, 'd_34_end': rr})
            else:
                submit_scan(f'{wd}/{str(step + 1)}_1/scan', charg, alpb,
                            kw={'concerted_trigger': False,
                                'at1': at_1, 'at2': at_2 + nat_1, 'd1': rr * 5, 'd2': rr, })
        else:
            print("gap_l wrong: ", reactants, push_dict_1)
    else:
        # full search
        candidates = kwargs.get('intercands')
        if not radical_trigger:
            gap_l = sorted(list(push_dict_1.keys()))[:candidates]
        else:
            gap_l = sorted(list(push_dict_1.keys()))[:int(candidates / 2)] + sorted(list(push_dict_2.keys()))[
                                                                             :int(candidates / 2)]
            push_dict_1.update(push_dict_2)
        count = 1
        if gap_l:
            for i in gap_l:
                # get reactants name and active atoms
                struc_1_nr, struc_2_nr = ast.literal_eval(push_dict_1[i][0])[0], ast.literal_eval(push_dict_1[i][0])[1]
                struc_1, struc_2 = reactants[int(struc_1_nr) - 1], reactants[int(struc_2_nr) - 1]
                at_l = ast.literal_eval(push_dict_1[i][1])
                at_1, at_2 = at_l[0][0], at_l[0][1]
                at_1, at_2 = int(at_1), int(at_2)

                coord_1 = rec[struc_1]['coord']
                nat_1 = coord_1.shape[0]
                coord_2 = rec[struc_2]['coord']
                nreact = struc_1.split('/')[0] + str(count)
                reactants_dict[nreact] = struc_1, struc_2
                coord_merge, rr = alignment(coord_1, coord_2, at_1, at_2)
                creat_folder(f'{wd}/{nreact}/scan')
                write_xyz(coord_merge, 'coord')
                charg = int(rec[struc_1]['charge'] + rec[struc_2]['charge'])

                if len(at_l) > 1:
                    at_3, at_4 = int(at_l[1][0]), int(at_l[1][1])
                    r3 = Coord(eval(f'coord_merge')).get_punkt_coord(at_3)
                    r4 = Coord(eval(f'coord_merge')).get_punkt_coord(at_4 + nat_1)
                    r = np.linalg.norm(r3 - r4)

                    submit_scan(f'{wd}/{nreact}/scan', charg, alpb,
                                kw={'at3': at_3, 'at4': at_4 + nat_1, 'concerted_trigger': True,
                                    'at1': at_1, 'at2': at_2 + nat_1, 'd1': rr * 5, 'd2': rr,
                                    'd_34_start': r, 'd_34_end': rr})

                else:
                    submit_scan(f'{wd}/{nreact}/scan', charg, alpb,
                                kw={'concerted_trigger': False,
                                    'at1': at_1, 'at2': at_2 + nat_1, 'd1': rr * 5, 'd2': rr, })

                count += 1
    return reactants_dict


def pull(reactants, rec, wd, step=0, **kwargs):
    kw = kwargs.get('kwargs')
    alpb = kw.get('alpb', None)
    reactor = kw.get('reactor_like', False)
    hydrogen = kw.get('hdrogen', True)

    reactants_dict = {}
    if not hydrogen:
        pull_dict = spr_intra_new(reactants, wd,hydrogen=False)
    pull_dict = spr_intra_new(reactants, wd)
    if pull_dict:

        # detemine active molecule and atoms
        min_gap = min(list(pull_dict.keys()))
        dis_struc = pull_dict[min_gap][0]
        dis_at_1 = pull_dict[min_gap][1]
        dis_at_2 = pull_dict[min_gap][2]
        coord = rec[dis_struc]['coord']
        C = Coord(coord)
        r1 = C.get_punkt_coord(dis_at_1)
        r2 = C.get_punkt_coord(dis_at_2)
        r = np.linalg.norm(r1 - r2)
        charg = rec[dis_struc]['charge']
        if reactor:
            # get reactant name and active atoms
            reactants_dict[f'{str(step + 1)}_2'] = dis_struc

            # write reactant xyz
            # creat_folder(f'{wd}/{str(step + 1)}_2')
            # write_xyz(coord, f'{dis_struc}')

            # submit scan job
            creat_folder(f'{wd}/{str(step + 1)}_2/scan')
            write_xyz(coord, 'coord')

            submit_scan(f'{wd}/{str(step + 1)}_2/scan', charg, alpb,
                        kw={'concerted_trigger': False,
                            'at1': dis_at_1, 'at2': dis_at_2, 'd1': r, 'd2': r * 5})

        else:
            nreact_count = str(kw.get('association_count', 5) + 1)
            nreact = dis_struc.split('/')[0] + nreact_count
            reactants_dict[nreact] = dis_struc
            creat_folder(f'{wd}/{nreact}/scan')
            write_xyz(coord, 'coord')

            submit_scan(f'{wd}/{nreact}/scan', charg, alpb,
                        kw={'concerted_trigger': False,
                            'at1': dis_at_1, 'at2': dis_at_2, 'd1': r, 'd2': r * 5})

    return reactants_dict
