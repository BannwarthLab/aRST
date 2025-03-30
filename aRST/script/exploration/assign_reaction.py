import random
import string
from collections import Counter
from os.path import join, exists
from aRST.script.toolbox import read_afo_csv, FragmentAnalyzer
from aRST.script.exploration.prepare_afo import getAFO
from aRST.script.callsys import CallORCA, CallMolBar
from aRST.script.exploration.interA import *
from aRST.script.exploration.intraA import *
from aRST.script.exploration.intraD import *
from aRST.script.head import print_simulation_details, print_jobhead
from aRST.script.exploration.network import rxnnetwork


class SingleCycleSearch:
    def __init__(self, reactants):
        self.reactants = reactants
        self.prepare_AFO()

    @staticmethod
    def check_featrues(struc):
        # check if coord, nel, afo, mbo exist
        coord = storage.allstruc.get_value(f"{struc}.struc_info.general.coord")
        nel = storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
        c = storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
        if not all(v is not None for v in [coord, nel, c]):
            print(f"Warning: {struc} doesn't have enough coord, nel, c, will be neglet in the reaction searching")

        afo = storage.allstruc.get_value(f"{struc}.struc_info.reactivity.afo")
        if afo is None:
            return False
        else:
            return True

    def prepare_AFO(self):
        for struc in self.reactants:
            if not self.check_featrues(struc):
                coord = storage.allstruc.get_value(f"{struc}.struc_info.general.coord")
                nel = storage.allstruc.get_value(f"{struc}.struc_info.general.nel")
                c = storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
                if coord.shape[0] == 1 and coord.iloc[0, 0] == "H":
                    mbo = None
                    afo = pd.DataFrame(
                        columns=['index', 'HOAO_a (eV)', 'HOAO_b (eV)', 'LUAO_a (eV)', 'LUAO_b (eV)', 'Atom'])
                    if c == 1:
                        # H+
                        afo.loc[0, :] = [0, -np.inf, -np.inf, -np.inf, -np.inf, 1]

                    elif c==-1:
                        # H-
                        afo.loc[0, :] = [0, np.inf, np.inf, np.inf, np.inf, 1]
                    else:
                        # H radical
                        afo.loc[0, :] = [0, np.inf, -np.inf, np.inf, -np.inf, 1]
                    storage.allstruc.set_value(f"{struc}.struc_info.reactivity.afo", afo)
                    storage.allstruc.set_value(f"{struc}.struc_info.reactivity.bo", mbo)


                else:
                    if str.lower(storage.setting.get_value("searching_setting.reactivity.afo_method")) == "gfn2":
                        afowd, status = getAFO(strucid=struc).get_from_xtb()
                    elif str.lower(storage.setting.get_value("searching_setting.reactivity.afo_method")) == "terachem":
                        afowd, status = getAFO(strucid=struc).get_from_tc()
                    else:
                        # defalut as orca sp
                        afowd, status = getAFO(strucid=struc).get_from_orca()
                    # read afo
                    afo = read_afo_csv(join(afowd, "ml_feature.csv"))

                    # ###### del H radical beta afo features (no beta electron at all) ######
                    # if nel == 1:
                    #     afo.loc[:, "HOAO_b (eV)"] = -np.inf
                    #     afo.loc[:, "LUAO_b (eV)"] = np.inf

                    storage.allstruc.set_value(f"{struc}.struc_info.reactivity.afo", afo)

                    # set bo
                    mbo = None
                    for file in [join(afowd, "wboe"), join(afowd, "bde"), join(afowd, "embo")]:
                        mbo = read_bo(file)
                        if mbo is not None:
                            break
                    storage.allstruc.set_value(f"{struc}.struc_info.reactivity.bo", mbo)

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
                if refreactants == at_l == refatl:
                    print(f"{rxntype} path guess for struc {structures}, "
                          f"atl {at_l} alreadt exist in rxn {rxnnetwork.get_node_index(node)}, jump")
                    return True
        return False

    def explore(self, strnr="0"):
        searchcount = 0
        output_interA, output_D, output_intraA = {}, {}, {}

        ## interA
        if storage.setting.get_value("searching_setting.inter", 1) and len(self.reactants) > 1:
            print()
            print_jobhead("Bimolecular Association")
            react_l = get_interA_react_l(self.reactants)
            scanmode = storage.setting.get_value("searching_setting.general.scan_mode", 0)
            ninter = storage.setting.get_value("searching_setting.inter.number", 1)
            if ninter == "all":
                ninter = len(react_l)
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
        if storage.setting.get_value("searching_setting.intra"):
            print()
            print_jobhead("Intramolecular Association")
            react_l = get_intraA_react_l(self.reactants)
            scanmode = storage.setting.get_value("searching_setting.general.scan_mode", 0)
            nintra = storage.setting.get_value("searching_setting.intra.number", 1)
            if nintra == "all":
                nintra = len(react_l)
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
                        res, searchcount = intrarmbff(self.reactants, entry, strnr, searchcount)
                        output_intraA.update(res)

        ## intraD
        if storage.setting.get_value("searching_setting.disso"):
            print()
            print_jobhead("Dissociation")
            react_l = get_intraD_react_l(self.reactants)
            scanmode = storage.setting.get_value("searching_setting.general.scan_mode", 0)
            ndisso = storage.setting.get_value("searching_setting.disso.number", 1)
            if ndisso == "all":
                ndisso = len(react_l)
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
        return out

    @staticmethod
    def follow_up_treatment_gfn(key, tmp_path_l):

        # doing sep opt for inter_fragments
        # output_d[count] = [{reactants info},{endxyz info}]
        reactants = tmp_path_l[0]["reactants"]
        cut = 4 if storage.setting.get_value("searching_setting.general.mb_level") == "loose" else 7
        reactants_mbl = [storage.allstruc.get_value(f"{i}.struc_info.general.mb")[:cut] for i in reactants]
        endxyz_info = tmp_path_l[1]

        # get constrained opt nfragment
        fragmentation = FragmentAnalyzer(endxyz_info["coord"], endxyz_info["bo"], endxyz_info["atomic_charge"])
        frag_data = fragmentation.fragment_data
        statue = fragmentation.status
        if not statue:
            raise ValueError(f"{key} fragmentation has error in endxyz, check!!!")

        # reset endxyz info
        endxyz_info_new = []

        def treate_more_geoms_in_endxyz(frag_data):
            # unconstrained opt for each fragment
            for frag in frag_data:
                coord = frag["coord"]
                charge = frag["charge_int"]
                nel = frag["nel"]
                atcharge = frag["atcharge"]

                random_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
                optwd = join(storage.setting.bufferwd, random_name, "xtb")
                create_folder(optwd)
                write_xyz(coord, "coord")
                if coord.shape[0] > 1:
                    uhf = 0 if nel % 2 == 0 else 1
                    status = CallxTB(optwd, custom_settings={"charge": charge, "uhf": uhf}).opt()
                    if status != 0:
                        status = CallxTB(optwd, custom_settings={"charge": charge, "uhf": uhf}).sp()
                    Esp_gfn = CallxTB(optwd).get_Esp_from_out()
                    coord = read_coord(join(optwd, "xtbopt.xyz"))
                    bo = read_bo(join(optwd, "wbo"))
                    constrain_l, ele_l = CallxTB(optwd).get_topo()
                    atcharge = read_atcharge(join(optwd, "charges"))

                    tested_fragemnation = FragmentAnalyzer(coord, bo, atcharge)
                    tested_frag_data = tested_fragemnation.fragment_data
                    if not tested_fragemnation.status:
                        raise ValueError(f"fragmentation has error, check {optwd}")
                    if len(tested_frag_data) > 2:
                        treate_more_geoms_in_endxyz(tested_frag_data)
                    else:

                        endxyz_info_new.append({"coord": coord,
                                                "atomic_charge": atcharge,
                                                "charge": charge,
                                                "nel": nel,
                                                "Esp_gfn": Esp_gfn,
                                                "strucwd": random_name,
                                                "atomic_connection": (constrain_l, ele_l), })
                else:
                    uhf = 0 if nel % 2 == 0 else 1
                    status = CallxTB(optwd, custom_settings={"charge": charge, "uhf": uhf}).sp()
                    Esp_gfn = CallxTB(optwd).get_Esp_from_out()
                    endxyz_info_new.append({"coord": coord,
                                            "atomic_charge": atcharge,
                                            "charge": charge,
                                            "nel": nel,
                                            "Esp_gfn": Esp_gfn,
                                            "strucwd": random_name,
                                            "atomic_connection": ([], [str(coord.iloc[0, 0])])})

        treate_more_geoms_in_endxyz(frag_data)

        # now evaluate if reaction simulation fall back to starting geom
        products_mbl = [CallMolBar(i["coord"]).get_mb() for i in endxyz_info_new]
        tmp_path_l.append(products_mbl)
        products_mbl = [_[:cut] for _ in products_mbl]
        bond_changed = True  # default as true
        if not "unknow" in products_mbl and not "unknow" in reactants_mbl:
            if Counter(products_mbl) == Counter(reactants_mbl):
                bond_changed = False
            else:
                bond_changed = True

        tmp_path_l[1] = endxyz_info_new
        tmp_path_l.append(bond_changed)

        # now output_d[count] = [{reactants info},[{P1}, {P2}],products_mbl, bond_changed_stuats]
        return tmp_path_l

    @staticmethod
    def follow_up_treatment_energy(tmp_path_l):
        ######## DFT OPT/SP ########
        dftopt = storage.setting.get_value("orca_setting.dftopt", False)

        for frag in tmp_path_l[1]:
            coord = frag["coord"]
            charge = frag["charge"]
            nel = frag["nel"]
            caltyp = "SP OPT" if coord.shape[0] > 1 and dftopt else "SP"
            m_l = [1, 3] if nel % 2 == 0 else [2]
            sp_dft, m_fin, coord_fin = None, None, None
            for m in m_l:
                wd = create_folder(
                    join(storage.setting.bufferwd, frag['strucwd'], "orca", str(m)))
                write_xyz(coord, "coord")
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
                        coord_fin = read_coord(join(wd, "orca.xyz")) if exists(join(wd, "orca.xyz")) else coord

            frag["Esp_orca"] = sp_dft
            frag["m"] = m_fin
            frag["coord"] = coord_fin

        ######## GFN TRV & SOLV ########
        for frag in tmp_path_l[1]:
            coord = frag["coord"]
            charge = frag["charge"]
            uhf = frag["m"] - 1
            etrv = 0
            if coord.shape[0] > 1:
                wd = create_folder(join(storage.setting.bufferwd, frag['strucwd'], "xtb", "trv"))
                write_xyz(coord, "coord")
                xtbprocess = CallxTB(wd, custom_settings={"charge": charge, "uhf": uhf, })
                staus = xtbprocess.trv()
                if staus == 0:
                    etrv = xtbprocess.get_Etrv()
            frag["Etrv_gfn"] = etrv

            esolv = 0
            if coord.shape[0] > 1:
                wd = create_folder(join(storage.setting.bufferwd, frag['strucwd'], "xtb", "solv"))
                write_xyz(coord, "coord")
                xtbprocess = CallxTB(wd, custom_settings={"charge": charge, "uhf": uhf, })
                staus = xtbprocess.solv()
                if staus == 0:
                    esolv = xtbprocess.get_Esolv()
            frag["Esolv_gfn"] = esolv
        return tmp_path_l

    def evalueate_path(self, key, tmp_path_l, layer=0):
        atl = tmp_path_l[0]["at_l"]
        simuhf = tmp_path_l[0]["simuhf"]
        # [{reactants info},[{P1}, {P2}],products_mbl, bond_changed_stuats]
        from collections import Counter
        reactants = tmp_path_l[0]["reactants"]
        unreactants = tmp_path_l[0]["unreacted"]
        simuhf = tmp_path_l[0]["simuhf"]

        cut = 4 if storage.setting.get_value("searching_setting.general.mb_level") == "loose" else 7
        reactants_mbl = [storage.allstruc.get_value(f"{i}.struc_info.general.mb")[:cut] for i in reactants]
        unreact_mbl = [storage.allstruc.get_value(f"{i}.struc_info.general.mb")[:cut] for i in unreactants]
        products_mbl = tmp_path_l[2]
        products_mbl = [_[:cut] for _ in products_mbl]

        info_temp = {"reactants": reactants,
                     "unreacted": unreactants,
                     "at_l": atl,
                     "sim_uhf": simuhf,
                     "products": [],
                     "scanwd": key,
                     # "E_rxn": E_rxn,
                     # "bond_changed": bondchange,
                     # "therdy_allow": therdy_allow,
                     "layer": layer,  # of which cycle step
                     }
        ###### determine the missing rxn info ######

        bondchange = tmp_path_l[3]
        if not bondchange:
            # reactants = products; adding new node with label bond not change
            products = reactants
            info_temp["products"] = products
            info_temp["bond_changed"] = False
            info_temp["E_rxn"] = 0
            info_temp["therdy_allow"] = True

            lastnode = rxnnetwork.get_lastnode_from_scanwd(key)
            edge_l = [lastnode]
            rxnnetwork.add_newnode(info=info_temp, edge_l=edge_l)
            print_simulation_details(status=1, key=key, atl=atl, sim_uhf=simuhf)

        else:
            info_temp["bond_changed"] = True
            # reactants != products; check if path guess thermodynamic assessable
            E_r = [storage.allstruc.get_value(f"{i}.struc_info.energy.sp_dft") for i in reactants]
            E_p = [i["Esp_orca"] for i in tmp_path_l[1]]
            if not any([_ is None for _ in E_r]) and not any([_ is None for _ in E_r]):
                E_rxn = sum(E_p) - sum(E_r)
            else:
                print(f"Warning: DFT energy calculation for scan {key} failed, change to gfn sp!")
                E_r = [storage.allstruc.get_value(f"{i}.struc_info.energy.sp_gfn") for i in reactants]
                E_p = [i["Esp_gfn"] for i in tmp_path_l[1]]
                E_rxn = sum(E_p) - sum(E_r)

            E_rxn_kcal = E_rxn * 627.503
            if E_rxn_kcal > storage.setting.get_value("searching_setting.general.elim", 300):
                info_temp["therdy_allow"] = False
                info_temp["E_rxn"] = E_rxn_kcal
                products = ["unknow"]
                info_temp["products"] = products  # not determing product mb
                rxnnetwork.add_newnode(info=info_temp)

            else:
                info_temp["therdy_allow"] = True
                info_temp["E_rxn"] = E_rxn_kcal
                # check if rxn exist already
                if not "unknow" in products_mbl and not "unknow" in reactants_mbl:
                    # compare mb
                    rxnexist = False
                    for node in rxnnetwork.graph.nodes:
                        info = rxnnetwork.get_reaction_info(node)
                        if int(info["layer"]) != 0:
                            reactants_ref = info["reactants"]
                            products_ref = info["products"]
                            unreact_ref = info["unreacted"]
                            refE_rxn = info["E_rxn"]
                            ref_rmbl = [storage.allstruc.get_value(f"{i}.struc_info.general.mb")[:cut] for i in
                                        reactants_ref]
                            ref_pmbl = [storage.allstruc.get_value(f"{i}.struc_info.general.mb")[:cut] for i in
                                        products_ref]
                            ref_umbl = [storage.allstruc.get_value(f"{i}.struc_info.general.mb")[:cut] for i in
                                        unreact_ref]
                            if Counter(products_mbl) == Counter(ref_pmbl) and \
                                    Counter(reactants_mbl) == Counter(ref_rmbl):
                                products = products_ref
                                # test unreactants too
                                if Counter(unreact_mbl) == Counter(ref_umbl):
                                    # not adding new node, but updating Erxn
                                    if refE_rxn is not None and E_rxn is not None and refE_rxn > E_rxn:
                                        rxnnetwork.graph.nodes[node]["info"]["E_rxn"] = E_rxn
                                        refnode_index = rxnnetwork.get_node_index(node)
                                        print_simulation_details(status=2, key=key, atl=atl,
                                                                 sim_uhf=simuhf, otherkey=refnode_index)
                                else:
                                    # adding as new node with new unreacted
                                    # add found new struc to storage
                                    products = products_ref
                                    info_temp["products"] = products
                                    rxnnetwork.add_newnode(info=info_temp)
                                    refnode_index = rxnnetwork.get_node_index(node)
                                    print_simulation_details(status=3, key=key, atl=atl, sim_uhf=simuhf,
                                                             otherkey=refnode_index)
                                rxnexist = True
                                break
                    if not rxnexist:
                        products = []
                        for i, pmb in enumerate(products_mbl):
                            add = True
                            if pmb != "unknow":
                                for struc in storage.allstruc.data.keys():
                                    refmb = storage.allstruc.get_value(f"{struc}.struc_info.general.mb")[:cut]
                                    if refmb == pmb:
                                        add = False
                                        products.append(struc)
                                        break
                            if add:
                                newstrucid = self.add_newstruc(tmp_path_l[1][i])
                                products.append(newstrucid)
                        info_temp["products"] = products
                        newnode_index = rxnnetwork.add_newnode(info=info_temp)
                        print_simulation_details(status=4, key=key, atl=atl, sim_uhf=simuhf, otherkey=newnode_index,
                                                 E_sp=E_rxn_kcal)
                else:
                    # can't identify, add as new node
                    products = []
                    for i, pmb in enumerate(products_mbl):
                        add = True
                        if pmb != "unknow":
                            for struc in storage.allstruc.data.keys():
                                refmb = storage.allstruc.get_value(f"{struc}.struc_info.general.mb")[:cut]
                                if refmb == pmb:
                                    add = False
                                    products.append(struc)
                                    break
                        if add:
                            newstrucid = self.add_newstruc(tmp_path_l[1][i])
                            products.append(newstrucid)
                    info_temp["products"] = products
                    newnode_index = rxnnetwork.add_newnode(info=info_temp)
                    print_simulation_details(status=6, key=key, atl=atl, sim_uhf=simuhf, otherkey=newnode_index,
                                             E_sp=E_rxn_kcal)

    @staticmethod
    def add_newstruc(strucinfo_d):
        newid = str(len(storage.allstruc.data.keys()))

        storage.allstruc.set_value(f"{newid}.struc_info.general.coord", strucinfo_d["coord"])
        storage.allstruc.set_value(f"{newid}.struc_info.general.charge", strucinfo_d["charge"])
        storage.allstruc.set_value(f"{newid}.struc_info.general.multiplicity", strucinfo_d["m"])
        storage.allstruc.set_value(f"{newid}.struc_info.general.nel", strucinfo_d["nel"])
        storage.allstruc.set_value(f"{newid}.struc_info.general.mb", CallMolBar(strucinfo_d["coord"]).get_mb())
        storage.allstruc.set_value(f"{newid}.struc_info.general.atcharge", strucinfo_d["atomic_charge"])
        storage.allstruc.set_value(f"{newid}.struc_info.general.atomic_connection", strucinfo_d["atomic_connection"])
        storage.allstruc.set_value(f"{newid}.struc_info.energy.sp_dft", strucinfo_d["Esp_orca"])
        storage.allstruc.set_value(f"{newid}.struc_info.energy.trv", strucinfo_d["Esp_gfn"])
        storage.allstruc.set_value(f"{newid}.struc_info.energy.solv", strucinfo_d["Esolv_gfn"])

        # move folder
        tmpwd = join(storage.setting.bufferwd, strucinfo_d["strucwd"])
        tarwd = join(storage.setting.strucwd, newid)
        if exists(tarwd):
            shutil.rmtree(tarwd)
        shutil.move(tmpwd, tarwd)
        return newid


def ReactionSearch(step=0):
    reactants_ini = storage.allstruc.reactants_ini
    refproducts_ini = storage.allstruc.refproducts
    print(f"Start doing reaction exploartions with given structues: ", reactants_ini)
    ncycle = storage.setting.get_value("searching_setting.general.cycle", 1)

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
            print_simulation_details(status=0)
            for itemd in explore_dl:
                for count, info_l in itemd.items():
                    # get_value products geom
                    info_l = singlecycle.follow_up_treatment_gfn(count, info_l)
                    # refine products info
                    info_l = singlecycle.follow_up_treatment_energy(info_l)
                    # now evaluate this reaction path
                    singlecycle.evalueate_path(count, info_l, layer=step + 1)
            print_simulation_details(status=5)
            # update reaction nexwork
            step += 1
        else:
            print("begin cycle ", step, ":")
            # cycle > 0
            nodes = rxnnetwork.get_allowed_layer_nodes(layer=step)
            nodes_scanwd = [rxnnetwork.get_reaction_info(_)["scanwd"] for _ in nodes]
            print("node = ", nodes_scanwd)
            for node,nodewd in zip(nodes,nodes_scanwd):
                print("starting exploring node ", nodewd)
                rxninfo = rxnnetwork.get_reaction_info(node)
                reactants_l = rxninfo["products"] + rxninfo["unreacted"]
                print("reactants for this node: ", reactants_l)
                prefix_rpid = rxninfo["scanwd"]
                singlecycle = SingleCycleSearch(reactants_l)
                explore_dl = singlecycle.explore(prefix_rpid)
                print_simulation_details(status=0)
                for itemd in explore_dl:
                    for count, info_l in itemd.items():
                        # get_value products geom
                        info_l = singlecycle.follow_up_treatment_gfn(count, info_l)
                        # refine products info
                        info_l = singlecycle.follow_up_treatment_energy(info_l)
                        # now evaluate this reaction path
                        # reaction_id = prefix_rpid +"_" + str(count.split("_")[-1])
                        singlecycle.evalueate_path(count, info_l, layer=step + 1)
                print_simulation_details(status=5)
            step += 1
