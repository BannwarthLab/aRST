import pandas as pd
import aRST.script.initialization as storage
from aRST.script.exploration.assign_reaction import SingleCycleSearch
from aRST.script.spr import SPR,SDP
from aRST.script.reactivity.print_stuff import *

def show_afo(printnumber=10):
    reactants = storage.allstruc.reactants_ini
    # preapre afo
    SingleCycleSearch(reactants)
    hydrogen = storage.setting.get_value("searching_setting.general.hydrogen", True)
    halogen = storage.setting.get_value("searching_setting.general.halogen", True)
    print(f"allow hydrogen atoms: {hydrogen}\nallow halogen atoms: {halogen}")
    spr = SPR(reactants, hydrogen, halogen)

    coord_l = [storage.allstruc.get_value(f"{_}.struc_info.general.coord") for _ in reactants]

    for i, list in enumerate([[spr.HOAO_a, spr.HOAO_b],[spr.LUAO_a, spr.LUAO_b]]):
        df_a,df_b =list
        df_a = df_a.apply(pd.to_numeric, errors='coerce')
        df_b = df_b.apply(pd.to_numeric, errors='coerce')
        if i ==0:
            type = "HOAO"
            df_a_sorted,df_b_sorted = df_a.stack().nlargest(printnumber).reset_index(), \
                                   df_b.stack().nlargest(printnumber).reset_index()

        else:
            type = "LUAO"
            df_a_sorted, df_b_sorted = df_a.stack().nsmallest(printnumber).reset_index(), \
                                       df_b.stack().nsmallest(printnumber).reset_index()
        printafo(status=0,printtype=type)

        at_l_a = df_a_sorted['level_0'].values
        struc_l_a = df_a_sorted['level_1'].values
        afo_l_a = df_a_sorted.values

        at_l_b = df_b_sorted['level_0'].values
        struc_l_b = df_b_sorted['level_1'].values
        afo_l_b = df_b_sorted.values

        for j in range(min(len(afo_l_a), len(afo_l_b))):
            struc_a = struc_l_a[j]
            at_a = at_l_a[j]
            ele_a =  coord_l[int(struc_a)].iloc[at_a,0]
            afo_a = afo_l_a[j][2]

            struc_b = struc_l_b[j]
            at_b = at_l_b[j]
            ele_b = coord_l[int(struc_b)].iloc[at_b, 0]
            afo_b = afo_l_b[j][2]


            printinfo = {"struc_a": struc_a, "at_a": at_a, "ele_a": ele_a, "afo_a": afo_a,
                         "struc_b": struc_b, "at_b": at_b, "ele_b": ele_b, "afo_b": afo_b}
            printafo(status=1,printtype=type, printinfo=printinfo)

        printafo(status=2)


    # # assess interA reactivity
    # if len(reactants)>1:
    #     print(f"Start doing reactivity assessment with given structues: ", reactants)
    #     hydrogen = storage.setting.get_value("searching_setting.general.hydrogen", True)
    #     halogen = storage.setting.get_value("searching_setting.general.halogen", True)
    #     print(f"allow hydrogen atoms: {hydrogen}; allow halogen atoms: {halogen}")
    #     spr = SPR(reactants,hydrogen,halogen)
    #
    #     for hadf,hbdf in zip(spr.HOAO_a,spr.HOAO_b):
    #         hadf_sorted, hbdf_sorted = hadf.stack().nlargest(printnumber).reset_index(),\
    #                                    hbdf.stack().nlargest(printnumber).reset_index()
    #
    #         HAat_index_l = hadf_sorted['level_0'].values
    #         HAstruc_l = hadf_sorted['level_1'].values
    #         HA_l = hadf_sorted.values
    #
    #         HBat_index_l = lbdf_sorted['level_0'].values
    #         HBstruc_l = hbdf_sorted['level_1'].values
    #         HB_l = hbdf_sorted.values
    #
    #         for i in range(min(len(HB_l),len(HA_l))):
    #             # HA
    #             struc_a = HAstruc_l[i]
    #             at_a = HAat_index_l[i]
    #             hoao_a = HA_l[i]
    #
    #             # HB
    #             struc_b = HBstruc_l[i]
    #             at_b = HBat_index_l[i]
    #             hoao_b = HB_l[i]



