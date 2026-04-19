from os.path import exists
from aRST.script.toolbox import read_afo_csv, FragmentAnalyzer, get_random_folder_name,\
    are_same_mb_list, are_same_struc_list
from aRST.script.exploration.prepare_afo import getAFO
from aRST.script.callsys import CallORCA, CallMolBar, CallgXTB
from aRST.script.exploration.interA import *
from aRST.script.exploration.intraA import *
from aRST.script.exploration.intraD import *
from aRST.script.head import print_simulation_details, print_jobhead
from aRST.script.exploration.network import rxnnetwork


class SingleCycleSearch:
    def __init__(self, reactants):
        self.reactants = reactants
        # add time
        from datetime import datetime
        start_time = datetime.now()
        self.prepare_AFO()
        end_time = datetime.now()
        running_time = end_time - start_time
        total_seconds = running_time.total_seconds()
        print(f"prepare AFO features takes {int(total_seconds // 3600)} hour "
              f"{int((total_seconds % 3600) // 60)} mins {total_seconds % 60:.2f} sec", flush=True)

    @staticmethod
    def check_featrues(struc):
        # check if coord, nel, afo, mbo exist
        s = storage.allstruc.get(struc)
        coordinate = s.struc_info.general.coordinate
        nel = s.struc_info.general.nel
        c = s.struc_info.general.charge
        if not all(v is not None for v in [coordinate, nel, c]):
            print(
                f"Warning: {struc} doesn't have enough coord, nel, c, this struc will be neglet in the reaction searching")

        afo = s.struc_info.reactivity.afo
        if afo is None:
            return False
        else:
            return True

    def prepare_AFO(self, afosubwd=None):
        for struc in self.reactants:
            if not self.check_featrues(struc):
                s = storage.allstruc.get(struc)
                coordinate = s.struc_info.general.coordinate
                ele_l = s.struc_info.general.ele_l
                nel = s.struc_info.general.nel
                c = s.struc_info.general.charge
                if len(ele_l) == 1 and ele_l[0] == "H":
                    mbo = None
                    afo = pd.DataFrame(
                        columns=['index', 'HOAO_a (eV)', 'HOAO_b (eV)', 'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom'])
                    if c == 1:
                        # H+
                        afo.loc[0, :] = [0, -np.inf, -np.inf, -np.inf, -np.inf, 1]

                    elif c == -1:
                        # H-
                        afo.loc[0, :] = [0, np.inf, np.inf, np.inf, np.inf, 1]
                    else:
                        # H radical
                        afo.loc[0, :] = [0, np.inf, -np.inf, np.inf, -np.inf, 1]
                    s.struc_info.reactivity.afo = afo
                    s.struc_info.reactivity.bo = mbo


                else:
                    if str.lower(storage.setting.searching_setting.reactivity.afo_method) == "gfn2":
                        afowd, status = getAFO(strucid=struc).get_from_xtb(afosubwd)
                    elif str.lower(storage.setting.searching_setting.reactivity.afo_method) == "terachem":
                        afowd, status = getAFO(strucid=struc).get_from_tc(afosubwd)
                    else:
                        # defalut as orca sp
                        afowd, status = getAFO(strucid=struc).get_from_orca(afosubwd)

                    # read afo
                    afo = read_afo_csv(join(afowd, "ml_feature.csv"))

                    # ###### del H radical beta afo features (no beta electron at all) ######
                    # if nel == 1:
                    #     afo.loc[:, "HOAO_b (eV)"] = -np.inf
                    #     afo.loc[:, "LUAO_b (eV)"] = np.inf

                    s.struc_info.reactivity.afo = afo

                    # set bo
                    mbo = None
                    for file in [join(afowd, "wboe"), join(afowd, "bde"), join(afowd, "embo")]:
                        mbo = read_bo(file)
                        if mbo is not None:
                            break
                    s.struc_info.reactivity.bo = mbo

    @staticmethod
    def rxn_guess_exist(entry, rxntype="interA"):
        # jump if guessed path already exist in rxn record
        if "struc_combi" in entry.keys():
            struc1, struc2 = entry["struc_combi"]
            structures = [struc1, struc2]
        else:
            structures = entry["struc"]
        at_l = entry["at_combi"]
        for node in rxnnetwork.graph.nodes:
            info = rxnnetwork.get_reaction_info(node)
            if int(info["layer"]) != 0:
                refreactants = info["reactants"]
                refatl = info["at_l"]
                if structures ==refreactants and at_l == refatl:
                    print(f"{rxntype} path guess for struc {structures}, "
                          f"atl {at_l} alreadt exist in rxn {rxnnetwork.get_node_index(node)}, jump")
                    return True
        return False

    def explore(self, strnr="0"):
        searchcount = 0
        output_interA, output_D, output_intraA = {}, {}, {}

        ## interA
        if len(self.reactants) > 1 and (
                storage.setting.searching_setting.inter.number > 0 or storage.setting.searching_setting.inter.atom_pair):
            print()
            print_jobhead("Bimolecular Association")
            from datetime import datetime
            start_time = datetime.now()
            react_l = get_interA_react_l(self.reactants)
            end_time = datetime.now()
            running_time = end_time - start_time
            total_seconds = running_time.total_seconds()
            print(f"prepare AFO gap ranking takes {int(total_seconds // 3600)} hour "
                  f"{int((total_seconds % 3600) // 60)} mins {total_seconds % 60:.2f} sec", flush=True)
            scanmode = storage.setting.searching_setting.general.scan_mode
            ninter = storage.setting.searching_setting.inter.number
            if ninter == "all":
                ninter = len(react_l)
            print(f"Paths for Ranking 0-{ninter - 1} will be simulated: ")
            if scanmode == 0:
                # call xtb scan
                for entry in react_l[:ninter]:
                    if not self.rxn_guess_exist(entry, "interA"):
                        res, searchcount = interxtbscan(self.reactants, entry, strnr, searchcount)
                        output_interA.update(res)
            else:
                # call mbff
                for entry in react_l[:ninter]:
                    if not self.rxn_guess_exist(entry, "interA"):
                        res, searchcount = intermbff(self.reactants, entry, strnr, searchcount)
                        output_interA.update(res)

        ## intraA
        if storage.setting.searching_setting.intra.number > 0 or \
                storage.setting.searching_setting.intra.atom_pair:
            print()
            print_jobhead("Intramolecular Association")
            # Testing if depotonation is energetically favorable
            if storage.setting.searching_setting.intra.depotonation:
                # Testing if depotonation is energetically favorable
                from aRST.script.exploration.depotonation import updating_reacntans
                updating_reacntans(self.reactants)
                self.prepare_AFO(afosubwd="depotonation")
            elif storage.setting.searching_setting.intra.potonation:
                # Testing if depotonation is energetically favorable
                from aRST.script.exploration.depotonation import updating_reacntans
                updating_reacntans(self.reactants, job="potonation")
                self.prepare_AFO(afosubwd="potonation")
            react_l = get_intraA_react_l(self.reactants)
            scanmode = storage.setting.searching_setting.general.scan_mode
            nintra = storage.setting.searching_setting.intra.number
            if nintra == "all":
                nintra = len(react_l)
            print(f"Paths for Ranking 0-{nintra - 1} will be simulated: ")
            if scanmode == 0:
                # call xtb scan
                for entry in react_l[:nintra]:
                    if not self.rxn_guess_exist(entry, "intrarA"):
                        res, searchcount = intraxtbscan(self.reactants, entry, strnr, searchcount)
                        output_intraA.update(res)
            else:
                # call mbff
                for entry in react_l[:nintra]:
                    if not self.rxn_guess_exist(entry, "intraA"):
                        res, searchcount = intrambff(self.reactants, entry, strnr, searchcount)
                        output_intraA.update(res)

        ## intraD
        if storage.setting.searching_setting.disso.number > 0 or \
                storage.setting.searching_setting.disso.atom_pair:
            print()
            print_jobhead("Dissociation")
            # Testing if depotonation is energetically favorable
            if storage.setting.searching_setting.disso.depotonation:
                # Testing if depotonation is energetically favorable
                from aRST.script.exploration.depotonation import updating_reacntans
                updating_reacntans(self.reactants)
                self.prepare_AFO(afosubwd="depotonation")
            elif storage.setting.searching_setting.disso.potonation:
                # Testing if depotonation is energetically favorable
                from aRST.script.exploration.depotonation import updating_reacntans
                updating_reacntans(self.reactants, job="potonation")
                self.prepare_AFO(afosubwd="potonation")
            react_l = get_intraD_react_l(self.reactants)
            scanmode = storage.setting.searching_setting.general.scan_mode
            ndisso = storage.setting.searching_setting.disso.number
            if ndisso == "all":
                ndisso = len(react_l)
            print(f"Paths for Ranking 0-{ndisso - 1} will be simulated: ")
            if scanmode == 0:
                # call xtb scan
                for entry in react_l[:ndisso]:
                    if not self.rxn_guess_exist(entry, "intraD"):
                        res, searchcount = dissoxtbscan(self.reactants, entry, strnr, searchcount)
                        output_D.update(res)

            else:
                # call mbff
                for entry in react_l[:ndisso]:
                    if not self.rxn_guess_exist(entry, "intraD"):
                        res, searchcount = dissombff(self.reactants, entry, strnr, searchcount)
                        output_D.update(res)

        out = {}
        for outd in [output_interA, output_D, output_intraA]:
            if outd:
                out.update(outd)
        out = [{key: value} for key, value in sorted(out.items(), key=lambda x: int(x[0].split('_')[-1]))]
        print()
        return out

    @staticmethod
    def follow_up_treatment_gfn(key, tmp_path_l):
        bufferwd = storage.setting.searching_setting.general.bufferwd

        # doing sep opt for inter_fragments
        # output_d[count] = [{reactants info},{endxyz info},{reaction info}]
        reactants = tmp_path_l[0]["reactants"]
        reactants_mbl = [storage.allstruc.get(i).struc_info.general.mb for i in reactants]
        endxyz_info = tmp_path_l[1]

        # get constrained opt nfragment
        fragmentation = FragmentAnalyzer(endxyz_info["coordinate"], endxyz_info["ele_l"], endxyz_info["bo"],
                                         endxyz_info["atomic_charge"])
        frag_data = fragmentation.fragment_data
        status = fragmentation.status
        if status != 0:
            raise ValueError(f"{key} fragmentation has error in endxyz, check!!!")

        # reset endxyz info
        endxyz_info_new = []

        def _trate_more_geoms_in_endxyz(frag_data):
            # target: doing unconstrained opt + (short md) for each fragments
            for frag in frag_data:
                coordinate = frag["coordinate"]
                ele_l = frag["ele_l"]
                charge = frag["charge_int"]
                nel = frag["nel"]
                atcharge = frag["atcharge"]
                if len(ele_l) > 1:
                    _handle_multi_at_geom(coordinate, ele_l, charge, nel)
                else:
                    _handle_single_at_geom(coordinate, ele_l, charge, nel, atcharge)

        def _handle_single_at_geom(coordinate, ele_l, charge, nel, atcharge):
            random_name = get_random_folder_name()
            spwd = join(bufferwd, key, random_name, "xtb")
            create_folder(spwd)
            write_xyz(coordinate, ele_l, "coord")
            uhf = 0 if nel % 2 == 0 else 1

            # doing xtb sp in random folder
            xtbparse = CallxTB(spwd, custom_settings={"charge": charge, "uhf": uhf})
            status = xtbparse.sp()
            Esp_gfn = xtbparse.get_Esp_from_out()
            endxyz_info_new.append({"coordinate": coordinate,
                                    "ele_l": ele_l,
                                    "atomic_charge": atcharge,
                                    "charge": charge,
                                    "nel": nel,
                                    "Esp_gfn": Esp_gfn,
                                    "strucwd": key + "/" + random_name,
                                    "atomic_connection": []})

        def _handle_multi_at_geom(coordinate, ele_l, charge, nel):
            # 1. fully relax
            random_name = get_random_folder_name()
            optwd = join(bufferwd, key, random_name, "xtb")
            create_folder(optwd)
            write_xyz(coordinate, ele_l, "coord")
            uhf = 0 if nel % 2 == 0 else 1

            # doing xtb opt in random folder
            xtbparse = CallxTB(optwd, custom_settings={"charge": charge, "uhf": uhf})
            status = xtbparse.opt()
            if status != 0:
                status = xtbparse.sp()  # switch to sp
            Esp_gfn = xtbparse.get_Esp_from_out()
            # receive opt res
            if exists(join(optwd, "xtbopt.xyz")):
                coordinate, ele_l = read_coord(join(optwd, "xtbopt.xyz"))
            bo = read_bo(join(optwd, "wbo"))
            constrain_l = get_topo(coordinate, ele_l, bo)
            atcharge = read_atcharge(join(optwd, "charges"))
            # now testing if more_geoms again
            tested_fragemnation = FragmentAnalyzer(coordinate, ele_l, bo, atcharge)
            tested_frag_data = tested_fragemnation.fragment_data
            if tested_fragemnation.status != 0:
                raise ValueError(f"fragmentation has error, check {optwd}")
            if len(tested_frag_data) > 1:
                _trate_more_geoms_in_endxyz(tested_frag_data)
            elif len(tested_frag_data) == 1:
                # add geom directly
                endxyz_info_new.append({"coordinate": coordinate,
                                        "ele_l": ele_l,
                                        "atomic_charge": atcharge,
                                        "charge": charge,
                                        "nel": nel,
                                        "Esp_gfn": Esp_gfn,
                                        "strucwd": key + "/" + random_name,
                                        "atomic_connection": constrain_l})

        _trate_more_geoms_in_endxyz(frag_data)
        # check up if number of nel; c; nat are correct
        totc_2 = sum(np.round(frag["charge"]) for frag in endxyz_info_new)
        totnel_2 = sum(int(frag["nel"]) for frag in endxyz_info_new)
        totnat_2 = sum(len(frag["ele_l"]) for frag in endxyz_info_new)
        totc = endxyz_info["charge"]
        totnel = endxyz_info["nel"]
        totnat = len(endxyz_info["ele_l"])
        all_equal = (totc == totc_2) and (totnel == totnel_2) and (totnat == totnat_2)
        if not all_equal:
            raise ValueError(f"Fragmentation in supra molecule in path {key} failed, check!")

        # now evaluate if reaction simulation fall back to starting geom
        products_mbl = []
        for item in endxyz_info_new:
            mb = CallMolBar(item["coordinate"], item["ele_l"]).get_mb()
            item["mb"] = mb
            products_mbl.append(mb)
        bond_change = True  # default as true
        if not "unknow" in products_mbl and not "unknow" in reactants_mbl:
            bond_change = False if are_same_mb_list(products_mbl,reactants_mbl,cut=True) else True

        tmp_path_l[1] = endxyz_info_new
        tmp_path_l[2]["bond_change"] = bond_change
        # now output_d[count] = [{reactants info},[{P1}, {P2}],{reaction info}]
        return tmp_path_l

    @staticmethod
    def follow_up_treatment_energy(tmp_path_l):
        bufferwd = storage.setting.searching_setting.general.bufferwd
        ######## DFT OPT/SP ########
        dftopt = storage.setting.orca_setting.dftopt
        for frag in tmp_path_l[1]:
            coordinate = frag["coordinate"]
            ele_l = frag["ele_l"]
            charge = frag["charge"]
            nel = frag["nel"]
            caltyp = "SP OPT" if len(ele_l) > 1 and dftopt else "SP"
            m_l = [1, 3] if nel % 2 == 0 else [2]
            if storage.setting.searching_setting.general.sp_method == "dft":
                sp_dft, m_fin, coordinate_fin = np.inf, m_l[0], coordinate  # initial setting in case DFT sp failed
                for m in m_l:
                    wd = create_folder(
                        join(bufferwd, frag['strucwd'], "orca", str(m)))
                    write_xyz(coordinate, ele_l, "coord")
                    orcaprocess = CallORCA(wd,
                                           custom_settings={"charge": charge,
                                                            "m": m,
                                                            "caltyp": caltyp})
                    staus = orcaprocess.call()
                    if staus == 0:
                        energy = orcaprocess.get_energy()
                        energy = np.inf if energy is None else energy  # jump when not converge
                        if sp_dft is None or energy < sp_dft:
                            sp_dft = energy
                            m_fin = m
                            coordinate_fin = coordinate
                            if exists(join(wd, "orca.xyz")):
                                coordinate_fin, ele_l = read_coord(join(wd, "orca.xyz"))

                frag["Esp"] = sp_dft
                frag["m"] = m_fin
                frag["coordinate"] = coordinate_fin
            else:
                sp_gxtb, m_fin, coordinate_fin = np.inf, m_l[0], coordinate  # initial setting in case DFT sp failed
                for m in m_l:
                    wd = create_folder(
                        join(bufferwd, frag['strucwd'], "gxtb", str(m)))
                    write_xyz(coordinate, ele_l, "coord")
                    gxtbprocess = CallgXTB(wd, custom_settings={"charge": charge, "uhf": m - 1})
                    with open(".CHRG", "w") as file:
                        file.write(f"{charge}\n")
                    with open(".UHF", "w") as file:
                        file.write(f"{m - 1}\n")
                    staus = gxtbprocess.sp()
                    if staus == 0:
                        energy = gxtbprocess.get_Esp_from_out()
                        energy = np.inf if energy is None else energy  # jump when not converge
                        if sp_gxtb is None or energy < sp_gxtb:
                            sp_gxtb = energy
                            m_fin = m
                            coordinate_fin = coordinate
                frag["Esp"] = sp_gxtb
                frag["m"] = m_fin
                frag["coordinate"] = coordinate_fin

        ######## GFN TRV & SOLV ######## # not use for now
        for frag in tmp_path_l[1]:
            coordinate = frag["coordinate"]
            ele_l = frag["ele_l"]
            charge = frag["charge"]
            uhf = frag["m"] - 1
            etrv = 0
            # if coord.shape[0] > 1:
            #     wd = create_folder(join(storage.setting.bufferwd, frag['strucwd'], "xtb", "trv"))
            #     write_xyz(coord, "coord")
            #     xtbprocess = CallxTB(wd, custom_settings={"charge": charge, "uhf": uhf, })
            #     staus = xtbprocess.trv()
            #     if staus == 0:
            #         etrv = xtbprocess.get_Etrv()
            frag["Etrv_gfn"] = etrv

            esolv = 0
            # if coord.shape[0] > 1:
            #     wd = create_folder(join(storage.setting.bufferwd, frag['strucwd'], "xtb", "solv"))
            #     write_xyz(coord, "coord")
            #     xtbprocess = CallxTB(wd, custom_settings={"charge": charge, "uhf": uhf, })
            #     staus = xtbprocess.solv()
            #     if staus == 0:
            #         esolv = xtbprocess.get_Esolv()
            frag["Esolv_gfn"] = esolv
        return tmp_path_l

    def evalueate_path(self, key, tmp_path_l, layer=0, linkoldnode=True, printdetail=True):
        atl = tmp_path_l[0]["at_l"]
        # [{reactants info},[{P1}, {P2}],products_mbl, bond_changed_stuats]
        from collections import Counter
        reactants = tmp_path_l[0]["reactants"]
        unreactants = tmp_path_l[0]["unreacted"]
        simuhf = tmp_path_l[0]["simuhf"]

        product_items = [(item["mb"],item["charge"],item["m"]) for item in tmp_path_l[1]]

        info_temp = {"reactants": reactants,
                     "unreacted": unreactants,
                     "at_l": atl,
                     "sim_uhf": simuhf,
                     "products": [],
                     "scanwd": key,
                     # "E_rxn": E_rxn,
                     # "bond_change": bond_change,
                     # "therdy_allow": therdy_allow,
                     "layer": layer,  # of which cycle step
                     "type": tmp_path_l[2]["type"]
                     # "rxnid" =
                     }
        ###### determine the missing rxn info ######

        bond_change = tmp_path_l[2]["bond_change"]
        if not bond_change:
            # reactants = products; adding new node with label bond not change
            products = reactants
            info_temp["products"] = products
            info_temp["bond_change"] = False
            info_temp["E_rxn"] = 0
            info_temp["therdy_allow"] = True
            info_temp["reaction_exist"] = 0  # not a new rxn

        else:
            info_temp["bond_change"] = True
            # reactants != products; check if path guess thermodynamic assessable
            E_p = [i["Esp"] for i in tmp_path_l[1]]
            E_r = [storage.allstruc.get(i).struc_info.energy.sp for i in reactants]
            all_energies = E_p + E_r
            if any(e is None or not np.isfinite(e) for e in all_energies):
                E_rxn = np.inf
            else:
                E_rxn = sum(E_p) - sum(E_r)

            E_rxn_kcal = E_rxn * 627.503 if not np.isinf(E_rxn) else np.inf
            # print(key,E_rxn_kcal)
            if E_rxn_kcal > storage.setting.searching_setting.general.elim:
                info_temp["therdy_allow"] = False
                info_temp["E_rxn"] = E_rxn_kcal
                products = ["unknow"]
                info_temp["products"] = products  #
                info_temp["reaction_exist"] = 2  # new reaction

            else:
                info_temp["therdy_allow"] = True
                info_temp["E_rxn"] = E_rxn_kcal
                # check if rxn exist already
                if not any("unknow" in mbcut for mbcut, _, _ in product_items):
                    info_temp["reaction_exist"] = 2  # preset new rxn
                    data = storage.allstruc.data
                    for node in [_ for _ in rxnnetwork.graph.nodes if
                                 int(rxnnetwork.get_reaction_info(_)["layer"]) != 0]:
                        refinfo = rxnnetwork.get_reaction_info(node)

                        # check reactant
                        reactants_ref = refinfo["reactants"]
                        if Counter(info_temp["reactants"]) != Counter(reactants_ref):
                            continue

                        # check product
                        product_ref = refinfo["products"]
                        ref_product_items = []
                        for sid in product_ref:
                            general = data[sid].struc_info.general
                            ref_product_items.append((general.mb,general.charge,general.multiplicity))

                        if are_same_struc_list(ref_product_items,product_items,cut=True):
                            # reaction exist
                            rxnexist = True
                            unreactants_product_ref = refinfo["unreacted"]
                            if Counter(unreactants) == Counter(unreactants_product_ref):
                                refnode_index = rxnnetwork.get_node_index_form_nodeid(node)
                                info_temp["reaction_exist"] = 0  # rxn exist
                                info_temp["ref_id"] = refnode_index
                            else:
                                refnode_index = rxnnetwork.get_node_index_form_nodeid(node)
                                info_temp["reaction_exist"] = 1  # rxn exist but there is new unreactants
                                info_temp["ref_id"] = refnode_index
                            break
                # check P structure list
                products = []
                for i, p_item in enumerate(product_items):
                    mb, c, m = p_item
                    if not "unknow" in mb:
                        sid_exist, refsid = storage.allstruc._struc_exist(mb, c, m)
                        if not sid_exist:
                            sid = self.add_newstruc(tmp_path_l[1][i])
                            products.append(sid)
                        else:
                            products.append(refsid)
                            refE = storage.allstruc.get(refsid).struc_info.energy.sp
                            if tmp_path_l[1][i].get("Esp") < refE:
                                self.update_current_structure_info(tmp_path_l[1][i], refsid)
                    else:
                        products.append("unknow")

                    info_temp["products"] = products
                    if "unknow" in products:
                        info_temp["reaction_exist"] = 3

        if not info_temp["bond_change"]:
            lastnode = rxnnetwork.get_lastnode_from_scanwd(key)
            edge_l = [lastnode]
            newnode_index = rxnnetwork.add_newnode(info=info_temp, edge_l=edge_l, linkoldnode=linkoldnode)
        else:
            newnode_index = rxnnetwork.add_newnode(info=info_temp, linkoldnode=linkoldnode)
        info_temp["rxnid"] = newnode_index
        # E_r_rxn will be updated when final print out
        return info_temp

    def add_newstruc(self, strucinfo_d):
        sid, s = storage.allstruc.creat_new_sid()
        m = strucinfo_d["m"]
        c = strucinfo_d["charge"]
        mb = strucinfo_d["mb"]
        mb_cut = " | ".join(mb.split(" | ")[:storage.setting.searching_setting.general.mbcut_length])
        storage.allstruc.unique_struc[(mb_cut, c, m)] = sid
        entry = s.struc_info.general
        entry.coordinate = strucinfo_d["coordinate"]
        entry.ele_l = strucinfo_d["ele_l"]
        entry.charge = strucinfo_d["charge"]
        entry.multiplicity = strucinfo_d["m"]
        entry.nel = strucinfo_d["nel"]
        entry.mb = mb
        entry.atcharge = strucinfo_d["atomic_charge"]
        entry.atomic_connection = strucinfo_d["atomic_connection"]
        s.struc_info.energy.sp = strucinfo_d["Esp"]

        # storage.allstruc.set_value(f"{newid}.struc_info.energy.trv", strucinfo_d["Etrv_gfn"])
        # storage.allstruc.set_value(f"{newid}.struc_info.energy.solv", strucinfo_d["Esolv_gfn"])

        # move folder
        tmpwd = join(storage.setting.searching_setting.general.bufferwd, strucinfo_d["strucwd"])
        tarwd = join(storage.setting.searching_setting.general.strucwd, str(sid))
        if exists(tarwd):
            shutil.rmtree(tarwd)
        shutil.move(tmpwd, tarwd)
        return sid

    @staticmethod
    def update_current_structure_info(strucinfo_d, refid):
        storage.allstruc.data[refid].struc_info.general.coordinate = strucinfo_d["coordinate"]
        storage.allstruc.data[refid].struc_info.general.atcharge = strucinfo_d["atomic_charge"]
        storage.allstruc.data[refid].struc_info.general.atomic_connection = strucinfo_d["atomic_connection"]
        storage.allstruc.data[refid].struc_info.energy.sp = strucinfo_d["Esp"]
        # move folder
        tmpwd = join(storage.setting.searching_setting.general.bufferwd, strucinfo_d["strucwd"])
        tarwd = join(storage.setting.searching_setting.general.strucwd, str(refid))
        if exists(tarwd):
            shutil.rmtree(tarwd)
        shutil.move(tmpwd, tarwd)


def ReactionSearch(step=0):
    reactants_ini = storage.allstruc.reactants_ini
    # refproducts_ini = storage.allstruc.refproducts
    print(f"Start doing reaction exploartions with given structues: ", reactants_ini)
    ncycle = storage.setting.searching_setting.general.cycle

    while step < ncycle:
        # find reaction pathes to be explored in this cycle
        if step == 0:
            print("begin cycle 0:")
            print()
            # cycle 0

            reactants_l = reactants_ini
            # add reactants_ini as 0 node
            rxnnetwork.add_newnode(info={"reactants": reactants_l,
                                         "unreacted": [],
                                         "products": [],
                                         "layer": 0})

            singlecycle = SingleCycleSearch(reactants_l)
            explore_dl = singlecycle.explore()
            print_info = []
            for itemd in explore_dl:
                for scan_wd, info_l in itemd.items():
                    # get_value products geom
                    info_l = singlecycle.follow_up_treatment_gfn(scan_wd, info_l)
                    # refine products info
                    info_l = singlecycle.follow_up_treatment_energy(info_l)
                    # now evaluate this reaction path
                    print_item = singlecycle.evalueate_path(scan_wd, info_l, layer=step + 1)
                    print_info.append(print_item)
            # before print out, update E_r_rxn
            print_info = rxnnetwork.self_check(print_info)
            print_simulation_details(print_info, "0")
            # update reaction nexwork
            step += 1
        else:
            print("begin cycle ", step, ":")
            # cycle > 0
            nodes = rxnnetwork.get_allowed_layer_nodes(layer=step)
            nodes_scanwd = [rxnnetwork.get_reaction_info(_)["scanwd"] for _ in nodes]
            print("node = ", nodes_scanwd)
            for node, nodewd in zip(nodes, nodes_scanwd):
                print("starting exploring node ", nodewd)
                rxninfo = rxnnetwork.get_reaction_info(node)
                reactants_l = rxninfo["products"] + rxninfo["unreacted"]
                print("reactants for this node: ", reactants_l)
                prefix_rpid = rxninfo["scanwd"]
                singlecycle = SingleCycleSearch(reactants_l)
                explore_dl = singlecycle.explore(prefix_rpid)
                print_info = []
                for itemd in explore_dl:
                    for count, info_l in itemd.items():
                        # get_value products geom
                        info_l = singlecycle.follow_up_treatment_gfn(count, info_l)
                        # refine products info
                        info_l = singlecycle.follow_up_treatment_energy(info_l)
                        # now evaluate this reaction path
                        # reaction_id = prefix_rpid +"_" + str(count.split("_")[-1])
                        print_item = singlecycle.evalueate_path(count, info_l, layer=step + 1)
                        print_info.append(print_item)
                # before print out, update E_r_rxn
                print_info = rxnnetwork.self_check(print_info)
                print_simulation_details(print_info, nodewd)
            step += 1
    # clean up buffer_struc
    if exists(storage.setting.searching_setting.general.bufferwd):
        import shutil
        shutil.rmtree(storage.setting.searching_setting.general.bufferwd, ignore_errors=True)
