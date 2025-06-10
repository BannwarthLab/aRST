import datetime


def print_head():
    program_name = "automated Reaction Search Tool (aRST) V2.0"
    timestamp = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    tot_width = 80
    border_width = 60
    side_padding = (tot_width - border_width) // 2
    side = " " * side_padding
    border = "=" * border_width

    print(side + border)
    print(side + f"|{' ' * (border_width - 2)}|")
    print(side + f"|{' ' * (border_width - 2)}|")
    print(side + f"|{program_name.center(border_width - 2)}|")
    print(side + f"|{' ' * (border_width - 2)}|")
    print(side + f"|{' ' * (border_width - 2)}|")
    print(side + border)

    print(side+"# Author:             Y. Chen")
    print(side+"# Institute:          Institute of Physical Chemistry")
    print(side+"#                     RWTH Aachen University")
    print(side+"# Version date:       03.2025")
    print(side+f"# Job running at:     {timestamp}")
    print(side + border)
    print()


def print_args(args):
    string = " "
    if args.explore:
        string = "Job type: reaction exploartion"
    elif args.afo_features:
        string = "Job type: reactivity assesment"
    elif args.align:
        string = "Job type: reactants alignment"
    print(string.center(80))


def print_parser(data):
    print("*read input:")
    for section, subkeys in data.items():
        # Start with section header
        line = f"[{section}]:"
        indent = " " * 4  # 4-space indent for continuation lines
        first_item_in_line = True  # Track if we're at the start of a new line

        # Add key-value pairs one by one, checking line length
        for k, v in subkeys.items():
            # Add comma only if it's not the first item in the line
            item = f", {k}={v}" if not first_item_in_line else f" {k}={v}"

            # If adding this item would exceed 60 chars, start a new line
            if len(line) + len(item) > 60:
                print(line)
                line = indent + f"{k}={v}"  # Start new line with indent (no comma)
                first_item_in_line = True  # Next item is the first in this line
            else:
                line += item
                first_item_in_line = False  # Next item is not the first

        print(line)  # Print the final line for the section
    print()
    print()


def print_simulation(reactants, at_l, charge_l, nel_l, reactionpath, count, totnel, totcharge, u, nat_1=None):
    totw = 80
    border_width = 60
    col1 = 15
    col2 = 15
    col3 = 15
    col4 = 15
    side_padding = (totw - border_width) // 2
    side = " " * side_padding

    print(side + "=" * border_width)
    title = f"Path Guess Scan {reactionpath}_{count}"
    print(title.center(totw))
    print(side + "-" * border_width)
    print(side + f"{'struc'.center(col1)}{'charge'.center(col2)}{'nel'.center(col3)}{'react atom'.center(col4)}")
    print(side + "-" * border_width)

    if len(at_l) > 1:
        for reactant, at1, at2, charge, nel in zip(reactants, at_l[0], at_l[1], charge_l, nel_l):
            atoms = f"({at1},{at2})"
            print(side + f"{str(reactant).center(col1)}"
                         f"{str(charge).center(col2)}"
                         f"{str(nel).center(col3)}"
                         f"{str(atoms).center(col4)}")

    else:
        for reactant, at, charge, nel in zip(reactants, at_l[0], charge_l, nel_l):
            print(side + f"{str(reactant).center(col1)}"
                         f"{str(charge).center(col2)}"
                         f"{str(nel).center(col3)}"
                         f"{str(at).center(col4)}")

    print(side + "-" * border_width)
    at1 = int(at_l[0][0])
    at2 = int(at_l[0][1])
    if nat_1:
        at2 += nat_1
    if len(at_l) > 1:
        at3 = int(at_l[1][0])
        at4 = int(at_l[1][1])
        if nat_1:
            at4 += nat_1
        string = f"({at1, at2},{at3, at4}"
    else:
        string = f"{at1, at2}"
    print(side + f"{f'P(uhf={u})'.center(col1)}"
                 f"{str(totcharge).center(col2)}"
                 f"{str(totnel).center(col3)}"
                 f"{str(string).center(col4)}")
    print(side + "=" * border_width)
    print()


def print_simulation_details(status,key=None, atl=None, sim_uhf=None, otherkey=None, E_sp=None):
    totw = 80
    border_width = 60
    col1 = 10
    col2 = 15
    col3 = 5
    col4 = 30
    side_padding = (totw - border_width) // 2
    side = " " * side_padding
    comments = None

    if status == 0:  # print head
        print(side + "=" * border_width)
        title = f"Reaction Simulation Summary"
        print(title.center(totw))
        print(side + "-" * border_width)
        print(side + f"{'scan':^{col1}} {'atl':^{col2}} {'uhf_sim':^{col3}} {'status':^{col4}}")
        print(side + "-" * border_width)
        return
    elif status == 5:  # print tail
        print(side + "=" * border_width)
        return
    elif status == 1:
        # reactants = products
        comments = "R=P"
    elif status == 2:
        # rxn exist
        comments = f"rxn exist in {otherkey}"

    elif status == 3:
        # rxn exist but different unreact
        comments = f" semi-new rxn: {otherkey}"
    elif status == 4:
        # new reaction
        comments = f" new rxn: {otherkey}, ΔE_sp ={E_sp:.2f}"
    elif status==6:
        # unknow in rmb or pmb
        comments = f" UMBnew rxn: {otherkey}, ΔE_sp ={E_sp:.2f}"
    print(side + f"{str(key).center(col1)}"
                 f"{str(atl).center(col2)}"
                 f"{str(sim_uhf).center(col3)}"
                 f"{str(comments).center(col4)}")


def print_jobhead(jobhead):
    string = "****"
    stringout = string + " " + jobhead + " " + string
    print(stringout.center(80))


def print_react_l(reactants,react_l,rtype="interA"):
    totw = 80
    border_width = 60
    col1 = 15
    col2 = 15
    col3 = 35
    side_padding = (totw - border_width) // 2
    side = " " * side_padding
    print(side + "*" * border_width)
    title = f"{rtype} Reactivity Assessment"
    print(title.center(totw))
    print()
    print(side + f"Reactants: {reactants}")
    print(side + f"tot.atom pairs: {len(react_l)}")
    print(side + "-" * border_width)
    print(side + f"{'struc pair':^{col1}} {'atom pair':^{col2}} {'AFO gap':^{col3}}")
    print(side + "-" * border_width)
    for entry in react_l[:20]:
        struc_l = entry["struc_combi"] \
            if "struc_combi" in entry.keys() else entry["struc"]
        at_l = entry["at_combi"]
        gap = entry["gap"]
        strgap = f"{gap:.2f}" if isinstance(gap, float) else "np.inf"
        print(side + f"{str(struc_l).center(col1)}"
                     f"{str(at_l).center(col2)}"
                     f"{strgap.center(col3)}")
    print(side + "*" * border_width)
