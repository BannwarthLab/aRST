import numpy as np

import aRST.script.initialization as storage
from aRST.script.exploration.assign_reaction import SingleCycleSearch
from aRST.script.exploration.network import rxnnetwork
from aRST.script.exploration.interA import *
from aRST.script.head import print_simulation_details,print_react_l
import random
from collections import Counter
from os.path import exists


def print_statistic(all_strucs_l):
    totw = 80
    border_width = 60
    side_padding = (totw - border_width) // 2
    side = " " * side_padding
    print(side + "*" * border_width)
    title = f"Current Summary"
    print(title.center(totw))
    print(side + "*" * border_width)
    print()
    print(side + f"tot. number of strucs = {len(all_strucs_l)}")
    print(side + f"tot. sorts of strucs = {len(set(all_strucs_l))}")
    print(side +f"top 10 structures with the highest number:")
    counter = Counter(all_strucs_l)
    top = counter.most_common(10)
    for item, count in top:
        print(side +f"struc {item}: {count}")
    print(side + "*" * border_width)
    print()

def print_selected_react_l(final_react_l):
    totw = 80
    border_width = 60
    col1 = 15
    col2 = 15
    col3 = 35
    side_padding = (totw - border_width) // 2
    side = " " * side_padding
    print(side + "*" * border_width)
    title = f"simulated reactions in this cycle"
    print(title.center(totw))
    print(side + "*" * border_width)
    for entry in final_react_l:
        struc_l = entry["struc_combi"] \
            if "struc_combi" in entry.keys() else entry["struc"]
        at_l = entry["at_combi"]
        gap = entry["gap"]
        strgap = f"{gap:.2f}" if isinstance(gap, float) else "np.inf"
        print(side + f"{str(struc_l).center(col1)}"
                     f"{str(at_l).center(col2)}"
                     f"{strgap.center(col3)}")
    print(side + "*" * border_width)
    print()


def select_react_l_interA(all_strucs_l, npair=10, gaplim=1, etemp=1000):
    reactants_l = list(set(all_strucs_l))
    # get full set of react_l
    react_l = get_interA_react_l(reactants_l,printall=False)
    # keep only item that afo gap < 1 eV
    screen_react_l = [item for item in react_l if item["gap"] < gaplim]
    if len(screen_react_l) <= npair:
        final_react_l = sorted(screen_react_l, key=lambda x: x["gap"])
    else:
        reactant_counts = Counter(all_strucs_l)
        nallstruc = len(all_strucs_l)
        weights = np.array([
            reactant_counts[item["struc_combi"][0]] * reactant_counts[item["struc_combi"][1]]/(nallstruc**2)
            for item in screen_react_l
        ])
        weights = weights / np.sum(weights)

        # use Metropolis-Hastings algo to sample npair of react_l
        count = 0
        k = 8.617e-5  # eV/K
        gap0 = screen_react_l[0]["gap"]  # init starting point with min gap
        final_react_l = [screen_react_l[0]]
        while count < npair-1:
            selected_index = np.random.choice(len(screen_react_l), p=weights)
            selected_pair = screen_react_l[selected_index]
            # selected_pair = random.choice(screen_react_l)
            if selected_pair in final_react_l:
                continue
            else:
                delta_gap = selected_pair["gap"] - gap0
                # np.exp(-E_a / (k_B * T))
                p = min(1, np.exp(-delta_gap / (k * etemp)))
                r = np.random.rand()
                if p >= r:
                    # accept
                    final_react_l.append(selected_pair)
                    gap0 = selected_pair["gap"]
                    count += 1
                else:
                    # reject
                    continue
        final_react_l = sorted(final_react_l, key=lambda x: x["gap"])
    # print out details
    print_selected_react_l(final_react_l)
    return final_react_l


def _get_path_possibility(newadded_node_l):
    # use bolzmann distribution for path possibility
    E_rxn_l = [rxnnetwork.get_reaction_info(rxnnode).get("E_rxn", 0) for rxnnode in newadded_node_l]
    k = 1.9872041 * 10 - 3  # kcal/mol
    etemp = storage.setting.get_value("xtb_setting.scan_etemp", 4000)
    p_l = np.array([np.exp(-(E / (k * etemp))) for E in E_rxn_l])
    print("!!")
    print(E_rxn_l)
    print(p_l)
    return p_l


def adjust_number_of_struc_in_reactor(newadded_node_l, layer, p_l):
    # self.number_of_struc = {"layer":[{"strucid": strucid,
    #                                   "num": num}, ...]}
    storage.allstruc.number_of_struc[str(layer)] = []
    for rxnnode in newadded_node_l:
        info = rxnnetwork.get_reaction_info(rxnnode)
        # check if this reaction accepted
        bond_changed = info["bond_changed"]
        if bond_changed:
            therdy_allow = info["therdy_allow"]
            if therdy_allow:
                reactants = info["reactants"]
                products = info["products"]
                p = p_l[newadded_node_l.index(rxnnode)]
                # check how many R struc can react
                current_num_of_rstruc = []
                for struc in reactants:
                    target_item = next(
                        (item for item in storage.allstruc.number_of_struc[str(layer)]
                         if item.get("strucid") == str(struc)), None)
                    if target_item:
                        current_num_of_rstruc.append(target_item["num"])
                    else:
                        oldterm = next(
                            (item for item in storage.allstruc.number_of_struc[str(int(layer) - 1)] if
                             item.get("strucid") == str(struc)), None)
                        current_num_of_rstruc.append(oldterm["num"])
                num = int(min(current_num_of_rstruc)*p)
                for struc in reactants:
                    target_item = next(
                        (item for item in storage.allstruc.number_of_struc[str(layer)]
                         if item.get("strucid") == str(struc)), None)
                    if target_item:
                        target_item["num"]-=num
                    else:
                        oldterm = next(
                            (item for item in storage.allstruc.number_of_struc[str(int(layer) - 1)] if
                             item.get("strucid") == str(struc)), None)
                        storage.allstruc.number_of_struc[str(layer)].append({"strucid":str(struc),"num":oldterm["num"]-num})
                for struc in products:
                    target_item = next(
                        (item for item in storage.allstruc.number_of_struc[str(layer)]
                         if item.get("strucid") == str(struc)), None)
                    if target_item:
                        target_item["num"] += num
                    else:
                        oldterm = next(
                            (item for item in storage.allstruc.number_of_struc[str(int(layer) - 1)] if
                             item.get("strucid") == str(struc)), None)
                        if oldterm:
                            storage.allstruc.number_of_struc[str(layer)].append({"strucid":str(struc),"num":oldterm["num"]+num})
                        else:
                            storage.allstruc.number_of_struc[str(layer)].append(
                                {"strucid": str(struc), "num":num})
    all_strucs_l = []
    for item in storage.allstruc.number_of_struc[str(layer)]:
        all_strucs_l.extend([item["strucid"]] * item["num"])
    return all_strucs_l

def pseudo_nanoreactor_inter(step=0):
    print("Start doing pseudo nanoreactor reaction explorations")
    all_strucs_l = storage.allstruc.reactants_ini
    ncycle = storage.setting.get_value("searching_setting.general.cycle", 1)
    ninter = storage.setting.get_value("searching_setting.inter.number", 1)
    react_etemp = storage.setting.get_value("xtb_setting.scan_etemp", 4000)
    scanmode = storage.setting.get_value("searching_setting.general.scan_mode", 0)
    while step < ncycle:
        print(f"begin cycle {step}:")
        print()
        output_interA = {}
        print_statistic(all_strucs_l)
        reactants_l = list(set(all_strucs_l))  # set of reactant structures

        # check if all reactant structure has AFO features
        afoprepare = SingleCycleSearch(reactants_l)
        reactl_interA = select_react_l_interA(reactants_l, ninter, gaplim=0.5, etemp=react_etemp)
        print_react_l(reactants_l, reactl_interA, rtype="interA")
        searchcount = 0

        if scanmode == 0:
            # call xtb scan
            for entry in reactl_interA:
                if not afoprepare.rxn_guess_exist(entry, "interA"):
                    reactants = entry["struc_combi"]
                    res, searchcount = interxtbscan(reactants, entry, str(step), searchcount)
                    output_interA.update(res)
        else:
            # call mbff
            for entry in reactl_interA:
                if not afoprepare.rxn_guess_exist(entry, "interA"):
                    reactants = entry["struc_combi"]
                    res, searchcount = intermbff(reactants, entry, str(step), searchcount)
                    output_interA.update(res)
        out = {}
        for outd in [output_interA]:
            if outd:
                out.update(outd)
        out = [{key: value} for key, value in sorted(out.items(), key=lambda x: int(x[0].split('_')[-1]))]
        print_simulation_details(status=0)
        old_node_l = list(rxnnetwork.graph.nodes)
        for itemd in out:
            for count, info_l in itemd.items():
                # get_value products geom
                info_l = afoprepare.follow_up_treatment_gfn(count, info_l)
                # refine products info
                info_l = afoprepare.follow_up_treatment_energy(info_l)
                # now evaluate this reaction path
                # reaction_id = prefix_rpid +"_" + str(count.split("_")[-1])
                afoprepare.evalueate_path(count, info_l, layer=step + 1,linkoldnode=False,printdetail=False)
        new_node_l = list(rxnnetwork.graph.nodes)
        newadded_node_l = [node for node in new_node_l if node not in set(old_node_l)]


        # update number of struc
        all_strucs_l = adjust_number_of_struc_in_reactor(newadded_node_l,layer=str(step+1),p_l=_get_path_possibility(newadded_node_l))
        print_simulation_details(status=5)
        step += 1

    if exists(storage.setting.bufferwd) and not (
    storage.setting.get_value("searching_setting.general.keepbuffstruc", False)):
        import shutil
        shutil.rmtree(storage.setting.bufferwd, ignore_errors=True)
