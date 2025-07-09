#!/usr/bin/python3.9
# import warnings
# warnings.simplefilter("ignore", UserWarning)
# warnings.simplefilter("ignore", RuntimeWarning)

import os
import shutil
import sys
import time
import traceback
from aRST.script.head import *
from argparse import ArgumentParser
from aRST.script.exploration.plot import plot_rxnnetwork
import aRST.script.initialization as init
from aRST.script.initialization import Setting




def main():

    parser = ArgumentParser(
        # name of program
        prog="aRST",
        # descriptionh
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
    parser.add_argument('--nanoreactor',
                        '-n',
                        action='store_true',
                        default=False,
                        help='Pseudo-Nanoreactor exploration based on given reactants.')
    parser.add_argument('--afo_features',
                        '-afo',
                        action='store_true',
                        default=False,
                        help='Print out reactive information of given structures.')
    parser.add_argument('--align',
                        '-a',
                        action='store_true',
                        default=False,
                        help='Aligning two strucs with given atom index.')
    parser.add_argument("-u", "--unbuffered", action="store_true",
                        help="Run in unbuffered mode (forces real-time output).")

    args = parser.parse_args()
    if args.unbuffered:
        sys.stdout = open(sys.stdout.fileno(), mode="w", buffering=1)
        sys.stderr = open(sys.stderr.fileno(), mode="w", buffering=1)

    print_head()

    # initialize
    setting = Setting(args.in_file)
    init.initialize(setting)

    print_parser(setting.config)
    print_args(args)

    if args.explore:
        from aRST.script.exploration.assign_reaction import ReactionSearch
        from datetime import datetime
        start_time = datetime.now()
        try:
            ReactionSearch()
            from aRST.script.toolbox import write_rxnout,write_strucout
            write_rxnout()
            write_strucout()
        except Exception as e:
            print(f"Error in main(): {e}", flush=True)
            traceback.print_exc()
        except SystemExit as e:
            print(f"SystemExit: {e}", flush=True)
            traceback.print_exc()
        end_time = datetime.now()
        running_time = end_time - start_time
        total_seconds = running_time.total_seconds()
        print()
        print(f"***** Job finished in {int(total_seconds // 3600)} hour "
              f"{int((total_seconds % 3600) // 60)} mins {total_seconds % 60:.2f} sec *****", flush=True)
        # plot_rxnnetwork()


    elif args.afo_features:
        from aRST.script.reactivity.afo_features import show_afo
        show_afo()
    elif args.align:
        from aRST.script.alignment.do_alignment import aligning_two_structures
        aligning_two_structures()
    elif args.nanoreactor:
        from aRST.script.nanoreactor.build_reaction import pseudo_nanoreactor_inter
        pseudo_nanoreactor_inter()




if __name__ == "__main__":

    main()

