import pandas as pd
from operation import creat_folder, write_xyz, xtb_ml
import time
import os


def ml(rec, reactants, pre_wd, xtbwd, alpb=None):
    if alpb is not None:
        alpb = [alpb][0]
    for name in reactants:
        if "ml_feature" not in rec[name].keys():
            coord = rec[name]['coord']
            creat_folder(f'{pre_wd}/{name}/ml')
            write_xyz(coord, 'coord')
            charg = rec[name]['charge']
            uhf = rec[name]['m']-1
            print(name, "doing ml  ")
            xtb_ml(xtbwd, charg, uhf, alpb)
            # double ml for H anion
            if coord.shape[0] == 1 and coord.iloc[0, 0] == 'H' and charg == -1:
                xtb_ml(xtbwd, charg,uhf, alpb)

    return rec


def collect_ml(rec, reactants, pre_wd):
    for name in reactants:
        if not name in rec.keys() or 'ml_feature' not in rec[name].keys():
            if os.path.isdir(f'{pre_wd}/{name}/ml'):
                if os.path.isfile(f'{pre_wd}/{name}/ml/xtb.out'):
                    ml_feature = pd.read_csv(f'{pre_wd}/{name}/ml/ml_feature.csv')
                    df = {}
                    if 'HOAO_a (eV)' in ml_feature.columns:
                        df['c'] = ml_feature['atomic partial charges']
                        df['cn'] = ml_feature['coordination number']
                        df['hoao_a'] = ml_feature['HOAO_a (eV)']
                        df['hoao_b'] = ml_feature['HOAO_b (eV)']
                        df['luao_a'] = ml_feature['LUAO_a (eV)']
                        df['luao_b'] = ml_feature['LUAO_b (eV)']
                        df['Atom'] = ml_feature['Atom']

                        ml_df = pd.DataFrame(df)
                    else:
                        df['c'] = ml_feature['atomic partial charges']
                        df['cn'] = ml_feature['coordination number']
                        df['hoao'] = ml_feature['HOAO (eV)']
                        df['luao'] = ml_feature['LUAO (eV)']
                        df['Atom'] = ml_feature['Atom']
                        ml_df = pd.DataFrame(df)
                    rec[name]['ml_feature'] = ml_df
        if not name in rec.keys() or 'wboe' not in rec[name].keys():
            try:
                if os.path.getsize(f'{pre_wd}/{name}/ml/wboe') != 0:
                    wboe = pd.read_csv(f'{pre_wd}/{name}/ml/wboe', header=None, sep='\s+')
                    rec[name]['wboe'] = wboe
                else:
                    rec[name]['wboe'] = []
            except:
                print(f'{pre_wd}/{name}/ml/wboe wrong')

    return rec
