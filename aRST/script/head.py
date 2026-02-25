import datetime


def print_head():
    program_name = "automated Reaction Search Tool (aRST) V2.2.0"
    timestamp = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    tot_width = 100
    border_width = 80
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

    print(side + "# Author:             Y. Chen, C. Bannwarth")
    print(side + "# Institute:          Institute of Physical Chemistry")
    print(side + "#                     RWTH Aachen University")
    print(side + "# Version date:       03.2026")
    print(side + f"# Job running at:     {timestamp}")
    print(side + border)
    print()


def print_args(args):
    string = " "
    if args.explore:
        string = "* Job type: Reaction Exploartion *"
    elif args.afo_features:
        string = "* Job type: Reactivity Assesment *"
    elif args.align:
        string = "* Job type: Reactants Alignment *"
    border_width = 80
    border = "=" * border_width
    print(border.center(100))
    print(string.center(100))
    print(border.center(100))
    print()


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


def print_simulation_details(info_l,node_name):
    totw = 100
    border_width = 90
    col1 = 18
    col2 = 10
    col3 = 12
    col4 = 18
    col5 = 10
    col6 = 32
    side_padding = (totw - border_width) // 2
    side = " " * side_padding

    # print head
    text_line = f"* Summary for the node {node_name} *"
    border = "*" * len(text_line)

    print(border.center(100))
    print(text_line.center(100))
    print(border.center(100))
    print()

    print(side + f"{'Path':^{col1}} {'Type':^{col2}} {'Reactant':^{col3}} {'Atom pair':^{col4}} "
                 f"{'UHF_sim':^{col5}} {'Status':^{col6}}")
    print(side + "-" * border_width)
    # info_temp = {"reactants": reactants,
    #                      "unreacted": unreactants,
    #                      "at_l": atl,
    #                      "sim_uhf": simuhf,
    #                      "products": [],
    #                      "scanwd": key,
    #                      # "E_rxn": E_rxn,
    #                      # "bond_change": bond_change,
    #                      # "therdy_allow": therdy_allow,
    #                      "layer": layer,  # of which cycle step
    #                      "type":tmp_path_l[2]["type"]
    #                      }
    for item in info_l:
        path= item["scanwd"]
        reactants=item["reactants"]
        formatted_reactants = ", ".join(f"[{','.join(map(str, inner))}]" for inner in item["reactants"])
        type=item["type"]
        formatted_atl = ", ".join(f"[{','.join(map(str, inner))}]" for inner in item["at_l"])
        uhf_sim = item["sim_uhf"]
        bond_change = item["bond_change"]
        therm_allow = item["therdy_allow"]
        E_sp = item["E_rxn"]
        reaction_exist_code = item["reaction_exist"]
        # 0: rxn exist; # 1: rxn exist but there is new unreactants
        # 2:new reaction; 3:unknow mb
        comments = ""
        if not bond_change:
            # reactants = products
            comments = "R=P"
        elif not therm_allow:
            # Erxn too large, add new node without determining products mb
            rxnid = item["rxnid"]
            comments = f"Unfavorable new rxn: {rxnid}, ΔE_sp ={E_sp:.2f} kcal/mol"
        else:
            rxnid = item["rxnid"]
            if reaction_exist_code==0:
                comments = f"rxn exist in {rxnid}"
            elif reaction_exist_code==1:
                comments = f"rxn exist in {rxnid}, but has different unreacted structures"
            elif reaction_exist_code == 2:
                comments = f"New rxn: {rxnid}, ΔE_sp ={E_sp:.2f}"
            elif reaction_exist_code == 3:
                comments = f"Unknown MB in product detected: {rxnid}, ΔE_sp ={E_sp:.2f}"
        print(side + f"{str(path).center(col1)}"
                     f"{str(type).center(col2)}"
                     f"{str(formatted_reactants).center(col3)}"
                     f"{str(formatted_atl).center(col4)}"
                     f"{str(uhf_sim).center(col5)}"
                     f"{str(comments).center(col6)}"
              )
    print(side + "-" * border_width)

def print_jobhead(jobhead):
    text_line = f"* {jobhead} *"
    border = "*" * len(text_line)

    print(border.center(100))
    print(text_line.center(100))
    print(border.center(100))
    print()


def print_react_l(reactants, react_l):
    totw = 100
    border_width = 80
    col1 = 10   # Ranking
    col2 = 15   # Struc combi
    col3 = 35   # Atom pair
    col4 = 15   # AFO gap / MBO / eMBO

    side_padding = (totw - border_width) // 2
    side = " " * side_padding

    title = "Reactivity Assessment"
    print(side + "-" * border_width)
    print(title.center(totw))
    print(side + f"Reactants: {reactants}, tot.atom pairs: {len(react_l)}")
    # print theo lvl

    print(side + "Top 10 reactivity predictions:")
    print(side + "-" * border_width)

    if "gap" in react_l[0]:
        gapkey, gaplabel = "gap", "AFO gap"
    elif "mbo" in react_l[0]:
        gapkey, gaplabel = "mbo", "MBO"
    else:
        gapkey, gaplabel = "embo", "eMBO"

    print(side + f"{'Ranking':^{col1}}{'Struc combi':^{col2}}{'Atom pair':^{col3}}{gaplabel:^{col4}}")
    print(side + "-" * border_width)

    for i, entry in enumerate(react_l[:10]):
        struc_l = entry.get("struc_combi", entry.get("struc", ""))
        at_l = entry["at_combi"]

        # at_l like [['5','0'],['2','4']] -> "[5,0], [2,4]"
        formatted_atl = ", ".join(f"[{','.join(map(str, inner))}]" for inner in at_l)

        gap = entry.get(gapkey)
        strgap = f"{gap:.2f}" if isinstance(gap, (int, float)) and gap != float("inf") else "np.inf"

        print(side +
              f"{i:^{col1}}"
              f"{str(struc_l):^{col2}}"
              f"{formatted_atl:^{col3}}"
              f"{strgap:^{col4}}")
    print(side + "-" * border_width)
    print(side +"The full AFO gap list is written to the file react_l.")
    print()
