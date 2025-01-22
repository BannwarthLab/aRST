#!/usr/bin/python3.9
# import warnings
# warnings.simplefilter("ignore", UserWarning)
# warnings.simplefilter("ignore", RuntimeWarning)

import os
import shutil
import time
import traceback

from argparse import ArgumentParser
from aRST.script.read import initial_struc_info, refine_struc_info,read_geom_d
from aRST.script.assign_reaction import ReactionSearch
from aRST.script.dictionaryInspector import DictionaryInspector


def main():

    parser = ArgumentParser(
        # name of program
        prog="aRST",
        # description
        description="an automated Reaction Search Tool",
        # bottom text
        epilog="By Ying"
    )

    # add arguments

    parser.add_argument("in_file", type=str, help="The input file.")
    parser.add_argument('--explore',
                        '-e',
                        action='store_true',
                        default=False,
                        help='Explore reaction pathes based on given reactants.')
    parser.add_argument('--reactor',
                        action='store_true',
                        default=False,
                        help='Nanoreactor likewise exploring thermpdynamic allowed products.')
    parser.add_argument('--ml_features',
                        '-ml',
                        action='store_true',
                        default=False,
                        help='Print out reactive information of given structures.')
    parser.add_argument('--align',
                        '-a',
                        action='store_true',
                        default=False,
                        help='Aligning two strucs with given atom index.')
    parser.add_argument('--opt',
                        action='store_true',
                        default=False,
                        help='Optimizing coord at gfn/dft level according to input file.')
    parser.add_argument('--time',
                        '-t',
                        action='store_true',
                        default=False,
                        help='Timing code runtime.')
    args = parser.parse_args()
    print(args)

    # read input
    import toml
    with open(args.in_file, "r") as f:
        setting_d = toml.load(f)
    print(setting_d)
    setting_inspertor = DictionaryInspector(setting_d)
    # determine the running folder
    if setting_inspertor.get_value("searching_setting.general.wd") is not None:
        wd = setting_inspertor.get_value("searching_setting.general.wd")
    else:
        wd = os.path.abspath(os.path.dirname(args.in_file))

    # read struc info and generate STRUC_rec
    reactants_ini, refproducts_ini, STRUC_rec, REFSTRUC_rec, REFSTRUC_rec,REACTION_rec=read_geom_d(wd,setting_d)
    STRUC_inspector = DictionaryInspector(STRUC_rec)
    REFSTRUC_inspector = DictionaryInspector(REFSTRUC_rec)


    if args.time:
        t1 = time.time()

    # determine task
    if args.explore:
        os.chdir(wd)
        # receive reaction path discovery record
        STRUC_rec, REACTION_rec = ReactionSearch(wd,reactants_ini,refproducts_ini,STRUC_rec, setting_d, {})
        #  examine if found ref reaction path
        from collections import Counter
        STRUC_inspector = DictionaryInspector(STRUC_rec)
        REF_inspector = DictionaryInspector(REFSTRUC_rec)
        ref_mb = [REF_inspector.get_value(f"{i}.struc_info.general.mb") for i in REFSTRUC_rec.keys()]

        triger = False
        for rxid in REACTION_rec.keys():
            unreacted_l = REACTION_rec[rxid][0]["unreacted"]
            products_l = REACTION_rec[rxid][0]["products"]
            products_l.extend(unreacted_l)
            pmb = [STRUC_inspector.get_value(f"{i}.struc_info.general.mb") for i in products_l]
            if Counter(pmb) == Counter(ref_mb):
                triger = True
                print("!!!! Congret, you have found ref. reaction path in reactionID:", rxid)
        if not triger:
            print("!!!! not found ref. reaction path with prod. mb", ref_mb)

        # write output file
        from aRST.script.dictionaryInspector import convert_df2json
        import json
        with open(f'{wd}/reaction.json', "w") as f:
            f.write(json.dumps(convert_df2json(REACTION_rec), indent=2))
            f.write('\n')

        with open(f'{wd}/struc.json', "w") as f:
            f.write(json.dumps(convert_df2json(STRUC_rec), indent=2))
            f.write('\n')
        if REFSTRUC_rec:
            with open(f'{wd}/refstruc.json', "w") as f:
                f.write(json.dumps(convert_df2json(REFSTRUC_rec), indent=2))
                f.write('\n')


    elif args.ml_features:
        from aRST.script.reactivity import get_reactivity
        os.chdir(wd)
        reactants_l = [i for i in STRUC_rec.keys()]
        if args.opt:
            # todo
            for strucid in STRUC_rec.keys():
                if STRUC_inspector.get_value(f"{strucid}.struc_info.general.sameas") is not None:
                    sameid = STRUC_inspector.get_value(f"{strucid}.struc_info.general.sameas")
                    STRUC_inspector.set_value(strucid, sameid)
                else:
                    STRUC_rec = refine_struc_info(wd, strucid, STRUC_rec, setting_d)
        # do stuff
        get_reactivity(wd, STRUC_rec, setting_d, reactants_l)
    elif args.reactor:
        os.chdir(wd)
        # refine STRUC_rec
        for strucid in STRUC_rec.keys():
            if STRUC_inspector.get_value(f"{strucid}.struc_info.general.sameas") is not None:
                sameid = STRUC_inspector.get_value(f"{strucid}.struc_info.general.sameas")
                STRUC_inspector.set_value(strucid, sameid)
            else:
                STRUC_rec = refine_struc_info(wd, strucid, STRUC_rec, setting_d)
        if REFSTRUC_rec:
            for strucid in REFSTRUC_rec.keys():
                if STRUC_inspector.get_value(f"{strucid}.struc_info.general.sameas") is not None:
                    sameid = STRUC_inspector.get_value(f"{strucid}.struc_info.general.sameas")
                    REFSTRUC_rec.set_value(strucid, sameid)
                else:
                    REFSTRUC_rec = refine_struc_info(wd, strucid, REFSTRUC_rec, setting_d)
        # receive reaction path discovery record
        try:
            STRUC_rec, REACTION_rec = ReactionSearch(wd, STRUC_rec, setting_d, {})

        except:
            traceback.print_exc()

        # write output file
        from aRST.script.dictionaryInspector import convert_df2json
        import json
        with open(f'{wd}/reaction.json', "w") as f:
            f.write(json.dumps(convert_df2json(REACTION_rec), indent=2))
            f.write('\n')

        with open(f'{wd}/struc.json', "w") as f:
            f.write(json.dumps(convert_df2json(STRUC_rec), indent=2))
            f.write('\n')


        pass
    elif args.align:
        from aRST.script.align import aligning
        os.chdir(wd)
        reactants_l = [i for i in STRUC_rec.keys()]
        if len(reactants_l) !=2:
            raise Exception("Number of input strucs is not 2, please check input")
        else:
            if args.opt:
                for strucid in STRUC_rec.keys():
                    if STRUC_inspector.get_value(f"{strucid}.struc_info.general.sameas") is not None:
                        sameid = STRUC_inspector.get_value(f"{strucid}.struc_info.general.sameas")
                        STRUC_inspector.set_value(strucid, sameid)
                    else:
                        STRUC_rec = refine_struc_info(wd, strucid, STRUC_rec, setting_d)
            aligning(wd,STRUC_rec,setting_d,reactants_l)

    if args.time:
        t2 = time.time()
        rt = t2-t1
        convert = time.strftime("%H:%M:%S", time.gmtime(rt))
        print(f"finished with {convert}")



if __name__ == "__main__":
    # import logging
    # log_to_stderr()
    # get_logger().setLevel(logging.DEBUG)
    main()
