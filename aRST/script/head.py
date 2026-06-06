import datetime
from typing import Any, Dict, List

import numpy as np

TOTAL_WIDTH = 100
DEFAULT_BORDER_WIDTH = 90


def _centered_banner(char: str, text: str, width: int = TOTAL_WIDTH) -> None:
    border = char * len(text)
    print(border.center(width))
    print(text.center(width))
    print(border.center(width))
    print()


def _framed_header(lines: List[str],
                   border_width: int = DEFAULT_BORDER_WIDTH,
                   width: int = TOTAL_WIDTH) -> None:
    side_padding = (width - border_width) // 2
    side = " " * side_padding
    border = "=" * border_width

    print(side + border)
    print(side + "|{:<{w}}|".format("", w=border_width - 2))
    print(side + "|{:<{w}}|".format("", w=border_width - 2))

    for line in lines:
        print(side + "|{:<{w}}|".format(line.center(border_width - 2), w=border_width - 2))

    print(side + "|{:<{w}}|".format("", w=border_width - 2))
    print(side + "|{:<{w}}|".format("", w=border_width - 2))
    print(side + border)


def print_head() -> None:
    program_name = "automated Reaction Search Tool (aRST) V2.2.0"
    timestamp = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    _framed_header(
        lines=[program_name],
        border_width=DEFAULT_BORDER_WIDTH,
        width=TOTAL_WIDTH,
    )

    side_padding = (TOTAL_WIDTH - DEFAULT_BORDER_WIDTH) // 2
    side = " " * side_padding
    border = "=" * DEFAULT_BORDER_WIDTH

    print(side + "# Author:             Y. Chen, C. Bannwarth")
    print(side + "# Institute:          Institute of Physical Chemistry")
    print(side + "#                     RWTH Aachen University")
    print(side + "# Version date:       03.2026")
    print(side + "# Job running at:     {}".format(timestamp))
    print(side + border)
    print()


def print_args(args: Any) -> None:
    job_type = ""
    if getattr(args, "explore", False):
        job_type = "* Job type: Reaction Exploration *"
    elif getattr(args, "afo_features", False):
        job_type = "* Job type: Reactivity Assessment *"
    elif getattr(args, "align", False):
        job_type = "* Job type: Reactants Alignment *"

    border_width = 80
    border = "=" * border_width
    print(border.center(TOTAL_WIDTH))
    if job_type:
        print(job_type.center(TOTAL_WIDTH))
    print(border.center(TOTAL_WIDTH))
    print()


def print_parser(data: Dict[str, Dict[str, Any]]) -> None:
    print("*read input:")
    for section, subkeys in data.items():
        line = "[{}]:".format(section)
        indent = " " * 4
        first_in_line = True

        for key, value in subkeys.items():
            item = ", {}={}".format(key, value) if not first_in_line else " {}={}".format(key, value)

            if len(line) + len(item) > 60:
                print(line)
                line = indent + "{}={}".format(key, value)
                first_in_line = True
            else:
                line += item
                first_in_line = False

        print(line)

    print()
    print()


def print_theo_setting() -> None:
    import aRST.script.initialization as storage

    setting = storage.setting

    width = 64  # inner box width (without outer centering)
    side_padding = (TOTAL_WIDTH - width) // 2
    side = " " * side_padding

    def box_line(text=""):
        print(side + "│ " + text.ljust(width - 4) + " │")

    def box_sep():
        print(side + "├" + "─" * (width - 2) + "┤")

    def box_top():
        print(side + "┌" + "─" * (width - 2) + "┐")

    def box_bottom():
        print(side + "└" + "─" * (width - 2) + "┘")

    def section(title):
        box_sep()
        box_line(title.center(width - 4))
        box_sep()

    def kv(key, value, key_width=18):
        box_line(f"{key:<{key_width}}: {value}")

    # Bind settings locally (cleaner)
    xtb = setting.xtb_setting
    orca = setting.orca_setting
    search = setting.searching_setting
    read = setting.read
    reactivity = search.reactivity

    # Prepare derived values
    dftmethod = ""
    dftsolvent = ""

    if search.general.sp_method == "dft":
        dftmethod = f"{orca.functional} {orca.basis}"
        dftsolvent = orca.cpcm

    scanmethod = "harmonic FF" if search.general.scan_mode == 1 else "xTB scan"

    # =============================
    # Start Box
    # =============================
    box_top()
    box_line("Theory Level Setting".center(width - 4))
    box_sep()

    # ---------------- Geometry ----------------
    section("Geometry Optimization")

    if not orca.dftopt:
        kv("Method", xtb.method)
        kv("Opt Level", xtb.optlevel or "normal")
        kv("Solvation", f"ALPB({xtb.alpb})" if xtb.alpb else "None")
        kv("Etemp", f"{xtb.struc_etemp} K")
    else:
        kv("Method", dftmethod)
        kv("Solvation", f"CPCM({dftsolvent})" if dftsolvent else "None")

    if read.general.jumprefine:
        kv("Re-opt input", "False")
    elif read.general.dftopt:
        kv("Re-opt input", dftmethod)
    elif read.general.gfn2opt:
        kv("Re-opt input", xtb.method)
    else:
        kv("Re-opt input", "False")

    # ---------------- Reactivity ----------------
    section("Reactivity Features")

    if reactivity.afo_method == "gfn2":
        kv("Method", "gfn2")
        kv("Solvation", f"ALPB({xtb.alpb})" if xtb.alpb else "None")

    elif reactivity.afo_method == "gxtb":
        kv("Method", "gxtb")

    elif reactivity.MO_method == "orca":
        afomethod = f"{reactivity.afo_orca_functional} {reactivity.afo_orca_basis}"
        kv("Method", afomethod)
        kv("Solvation", f"CPCM({reactivity.afo_orca_cpcm})"
        if reactivity.afo_orca_cpcm else "None")

    elif reactivity.MO_method == "terachem":
        afomethod = f"{reactivity.afo_tc_functional} {reactivity.afo_tc_basis}"
        kv("Method", afomethod)
        kv("Solvation", f"eps={reactivity.afo_tc_epsilon}"
        if reactivity.afo_tc_epsilon else "None")

    kv("AFO Gap tolerance", f"{search.general.equivalent_afogap_tol} eV")

    if search.general.exclude_elements:
        kv("Exclude elements", search.general.exclude_elements)

    if search.general.exclude_pairs:
        kv("Exclude pairs", search.general.exclude_pairs)

    # ---------------- Simulation ----------------
    section("Reaction Simulation")

    kv("Method", scanmethod)
    kv("Etemp", f"{xtb.scan_etemp} K")

    # ---------------- SP ----------------
    section("Single-Point Energy")

    kv("Method", search.general.sp_method)

    if search.general.sp_method == "dft":
        kv("Solvation", f"CPCM({dftsolvent})" if dftsolvent else "None")

    box_bottom()
    print()


def print_jobhead(jobhead: str) -> None:
    _centered_banner("*", "* {} *".format(jobhead), width=TOTAL_WIDTH)


def print_simulation_details(info_l: List[Dict[str, Any]], node_name: str, evalmode="g", tseval=False) -> None:
    _centered_banner("*", "* Exploration Summary of Node {} *".format(node_name), width=TOTAL_WIDTH)

    col_i = 4
    col_rxn = 6
    col_sw = 10
    col_r = 10
    col_ap = 14
    col_p = 12
    col_ts = 10
    col_erxn = 10
    col_note = 20

    side_padding = (TOTAL_WIDTH - 96) // 2
    side = " " * side_padding

    header = (
            side
            + "{:^{w1}} {:^{w2}} {:^{w3}} {:^{w4}} {:^{w5}} {:^{w6}} {:^{w7}} {:^{w8}} {:^{w9}}".format(
        "#", "RXNID", "Scan_wd", "R", "Atom Pair", "P", "ΔE‡", "ΔE", "Note",
        w1=col_i, w2=col_rxn, w3=col_sw, w4=col_r, w5=col_ap, w6=col_p, w7=col_ts, w8=col_erxn, w9=col_note
    )
    )
    if evalmode != "g":
        msg = "ΔE = Esp_gas(P)-Esp_gas(R)"
        if tseval:
            msg += ", ΔE‡ = Esp_gas(TS)-Esp_gas(R)"
    else:
        msg = "ΔE = G(P)-G(R)"
        if tseval:
            msg += ", ΔE‡ = G(TS)-G(R)"
        msg += ", G(i) = Esp_gas(i)+Gtrv(i)+Gsolv_corr(i)"

    msg += ", unit given in kcal/mol."

    print(side + msg)
    print(side + "-" * 96)
    print(header)
    print(side + "-" * 96)

    energy_key = "E_rxn" if evalmode != "g" else "G_rxn"
    TS_key = "dE_ts" if evalmode != "g" else "dG_ts"

    def safe_float(val):
        return float(val) if isinstance(val, (int, float)) else val
    for i, item in enumerate(info_l):
        reactants = item.get("reactants")
        products = item.get("products")
        scanwd = item.get("scanwd")
        atom_pairs = item.get("at_l")
        bond_change = bool(item.get("bond_change", False))
        therm_allow = bool(item.get("therdy_allow", True))
        rxnid = item.get("rxnid", "-")
        ref_id = item.get("ref_id", "-")
        reaction_exist = bool(item.get("reaction_exist", False))

        E_rxn = item.get(energy_key, "-")
        Ets_rxn = item.get(TS_key, "-")


        E_rxn = f"{E_rxn:.1f}" if isinstance(E_rxn, (int, float)) else str(E_rxn)
        Ets_rxn = f"{Ets_rxn:.1f}" if isinstance(Ets_rxn, (int, float)) else str(Ets_rxn)
        comment = item.get("note", "")

        line = (
            "{:^{w1}} {:^{w2}} {:^{w3}} {:^{w4}} {:^{w5}} {:^{w6}} {:^{w7}} {:^{w8}} {:<{w9}}".format(
                i, str(rxnid), str(scanwd), str(reactants), str(atom_pairs),
                str(products), safe_float(Ets_rxn),safe_float(E_rxn),
                str(comment),
                w1=col_i, w2=col_rxn, w3=col_sw, w4=col_r, w5=col_ap, w6=col_p, w7=col_ts, w8=col_erxn, w9=col_note
            )
        )
        print(side + line + side)

    print(side + "-" * DEFAULT_BORDER_WIDTH)

def print_population_summary(
    population: Dict[int, float],
    abundance: Dict[int, float],
    removed_info: List[Dict[str, Any]],
    step: int,
    tau: float,
    top_n: int = 10,
) -> None:
    _centered_banner(
        "*",
        "* Population / Abundance Summary of Layer {} *".format(step),
        width=TOTAL_WIDTH,
    )

    col_i = 4
    col_sid = 10
    col_abd = 14
    col_pop = 14
    col_reason = 40

    survivor_table_width = col_i + col_sid + col_abd + col_pop + 4
    removed_table_width = col_i + col_sid + col_abd + col_pop + col_reason + 5
    table_width = max(survivor_table_width, removed_table_width)

    side_padding = max((TOTAL_WIDTH - table_width) // 2, 0)
    side = " " * side_padding

    print(side + f"tau = {tau:.6e}")
    print(side + f"survivor number = {len(abundance)}")
    print(side + f"removed species number = {len(removed_info)}")
    print(side + "-" * table_width)

    # ============================================================
    # Surviving species, only show top 10
    # ============================================================
    print(side + f"Top {min(top_n, len(abundance))} surviving species:")
    print(side + "-" * survivor_table_width)

    survivor_header = "{:^{w1}} {:^{w2}} {:^{w3}} {:^{w4}}".format(
        "#",
        "SpeciesID",
        "Abundance",
        "Population",
        w1=col_i,
        w2=col_sid,
        w3=col_abd,
        w4=col_pop,
    )

    print(side + survivor_header)
    print(side + "-" * survivor_table_width)

    top_survivors = sorted(
        abundance.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    for i, (sid, abd) in enumerate(top_survivors[:top_n]):
        pop = population.get(sid, 0.0)

        line = "{:^{w1}} {:^{w2}} {:^{w3}} {:^{w4}}".format(
            i,
            str(sid),
            f"{abd:.1%}",
            f"{pop:.4e}",
            w1=col_i,
            w2=col_sid,
            w3=col_abd,
            w4=col_pop,
        )

        print(side + line)

    print(side + "-" * survivor_table_width)

    # Complete surviving species list
    survivor_ids = [
        sid for sid, _ in top_survivors
    ]

    print(side + "Complete surviving species list:")
    print(side + str(survivor_ids))
    print(side + "-" * table_width)

    # ============================================================
    # Removed species
    # ============================================================
    print(side + "Removed species:")
    print(side + "-" * removed_table_width)

    removed_header = "{:^{w1}} {:^{w2}} {:^{w3}} {:^{w4}} {:<{w5}}".format(
        "#",
        "SpeciesID",
        "Abundance",
        "Population",
        "Reason",
        w1=col_i,
        w2=col_sid,
        w3=col_abd,
        w4=col_pop,
        w5=col_reason,
    )

    print(side + removed_header)
    print(side + "-" * removed_table_width)

    if not removed_info:
        print(side + "No species removed.")
    else:
        for i, item in enumerate(removed_info):
            sid = item.get("sid", item.get("species", item.get("nodeid", "-")))

            abd = item.get("abundance", None)
            pop = item.get("population", population.get(sid, 0.0))
            reason = item.get("reason", item.get("note", "-"))

            if isinstance(abd, (int, float)):
                abd_str = f"{abd:.1%}"
            else:
                abd_str = "-"

            if isinstance(pop, (int, float)):
                pop_str = f"{pop:.4e}"
            else:
                pop_str = str(pop)

            line = "{:^{w1}} {:^{w2}} {:^{w3}} {:^{w4}} {:<{w5}}".format(
                i,
                str(sid),
                abd_str,
                pop_str,
                str(reason),
                w1=col_i,
                w2=col_sid,
                w3=col_abd,
                w4=col_pop,
                w5=col_reason,
            )

            print(side + line)

    print(side + "-" * table_width)


def print_react_l(reactants: Any, react_l: List[Dict[str, Any]]) -> None:
    totw = TOTAL_WIDTH
    border_width = 80
    col1 = 10  # Ranking
    col2 = 15  # Struc combi
    col3 = 35  # Atom pair
    col4 = 15  # AFO gap / MBO / eMBO

    side_padding = (totw - border_width) // 2
    side = " " * side_padding

    title = "Reactivity Assessment"
    print(side + "-" * border_width)
    print(title.center(totw))
    print(side + "Reactants: {}, tot.atom pairs: {}".format(reactants, len(react_l)))
    print(side + "Top 10 reactivity predictions:")
    print(side + "-" * border_width)

    gapkey = "embo"
    gaplabel = "eMBO"
    if react_l:
        if "gap" in react_l[0]:
            gapkey, gaplabel = "gap", "AFO gap"
        elif "mbo" in react_l[0]:
            gapkey, gaplabel = "mbo", "MBO"

    print(side + "{:^{c1}}{:^{c2}}{:^{c3}}{:^{c4}}".format(
        "Ranking", "Struc combi", "Atom pair", gaplabel,
        c1=col1, c2=col2, c3=col3, c4=col4
    ))
    print(side + "-" * border_width)

    for i, entry in enumerate(react_l[:10]):
        struc_l = entry.get("reactants", entry.get("reactants", ""))
        at_l = entry.get("atl", [])

        formatted_atl = ", ".join(
            "[{}]".format(",".join(map(str, inner))) for inner in at_l
        )

        gap = entry.get(gapkey)
        if isinstance(gap, (int, float)) and gap != float("inf"):
            strgap = "{:.2f}".format(gap)
        else:
            strgap = "np.inf"

        print(side + "{:^{c1}}{:^{c2}}{:^{c3}}{:^{c4}}".format(
            i, str(struc_l), formatted_atl, strgap,
            c1=col1, c2=col2, c3=col3, c4=col4
        ))

    print(side + "-" * border_width)
    print(side + "The full AFO gap list is written to the file react_l.")
    print()
