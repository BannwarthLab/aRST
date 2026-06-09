#!/usr/bin/python3.9
import sys
import traceback
from aRST.script.head import *
from argparse import ArgumentParser
import aRST.script.initialization as init
from aRST.script.initialization import Setting
from datetime import datetime




def main():

    parser = ArgumentParser(
        # name of program
        prog="aRST",
        # descriptionh
        description="an automated Reaction Search Tool",
        # bottom text
        epilog=""
    )

    # add arguments

    parser.add_argument("in_file", type=str,nargs="?", help="The input file.")
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
    parser.add_argument("-d", "--debug", action="store_true",
                        help="keep buffer folder")

    parser.add_argument(
        "--view",
        "-v",
        type=str,
        default=None,nargs="+",
        metavar="FILE, Strucjson",
        help="put in the to be viewed file and struc.json"
    )

    args = parser.parse_args()
    if args.unbuffered:
        sys.stdout = open(sys.stdout.fileno(), mode="w", buffering=1)
        sys.stderr = open(sys.stderr.fileno(), mode="w", buffering=1)

    if args.view is not None:
        if len(args.view) ==3:
            from aRST.script.visual.aseView import ViewSimPath
            ViewSimPath(args.view).show()
            return
        if len(args.view)==2 and not args.view[0].endswith("log"):
            from aRST.script.visual.aseView import ViewReactiveAtl
            ViewReactiveAtl(args.view).show()
            return



    print_head()

    # initialize
    setting = Setting(args.in_file)
    init.initialize(setting)

    # print_parser(setting.config)
    # print_args(args)

    if args.explore:
        from aRST.script.exploration.reactionsearch import ReactionSearch

        print_theo_setting()
        start_time = datetime.now()
        try:
            # prepare the input reactant
            reactants_ini = init.register_init_input_reactants(init.setting, init.allstruc)
            ReactionSearch(reactants_ini,debugtrigger=args.debug)

        except Exception as e:
            print(f"Error in main(): {e}", flush=True)
            traceback.print_exc()
        except SystemExit as e:
            print(f"SystemExit: {e}", flush=True)
            traceback.print_exc()
        finally:
            from aRST.script.toolbox import write_rxnout, write_strucout
            write_rxnout()
            write_strucout()
        end_time = datetime.now()
        running_time = end_time - start_time
        total_seconds = running_time.total_seconds()
        print()
        print(f"***** Job finished in {int(total_seconds // 3600)} hour "
              f"{int((total_seconds % 3600) // 60)} mins {total_seconds % 60:.2f} sec *****", flush=True)
        # clean up
        # from aRST.script.exploration.reactionsearch import clean_up_bufferwd
        # clean_up_bufferwd(args.debug, init.setting.searching_setting.general.bufferwd)
        # plot_rxnnetwork()


    elif args.afo_features:
        from aRST.script.reactivity.afo_features import show_afo
        show_afo()
    elif args.align:
        from aRST.script.alignment.do_alignment import aligning_two_structures
        reactants_ini = init.register_init_input_reactants(init.setting, init.allstruc)
        aligning_two_structures(reactants=reactants_ini)
    elif args.nanoreactor:
        print_theo_setting()
        start_time = datetime.now()
        try:
            from aRST.script.nanoreactor.worflow import Cancer
            reactants_ini = init.register_init_input_reactants(init.setting, init.allstruc)
            Cancer(reactants_ini,debugtrigger=args.debug)
        except Exception as e:
            print(f"Error in main(): {e}", flush=True)
            traceback.print_exc()
        except SystemExit as e:
            print(f"SystemExit: {e}", flush=True)
            traceback.print_exc()
        finally:
            from aRST.script.toolbox import write_rxnout, write_strucout
            write_rxnout()
            write_strucout()
        end_time = datetime.now()
        running_time = end_time - start_time
        total_seconds = running_time.total_seconds()
        print()
        print(f"***** Job finished in {int(total_seconds // 3600)} hour "
              f"{int((total_seconds % 3600) // 60)} mins {total_seconds % 60:.2f} sec *****", flush=True)




if __name__ == "__main__":

    main()

