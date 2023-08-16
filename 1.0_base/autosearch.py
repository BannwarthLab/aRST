import shutil
import os
import traceback
import function as F

# r_l = [['h', 'n2o'],['h', 'ch3f'], ['h', 'f2'], ['ch3', 'clf'], ['f-', 'ch3cl'], ['fch3clcomp1', 'fch3clcomp2'],
#      ['oh-', 'ch3f'], ['hoch3fcomp2', 'hoch3fcomp1'], ['h', 'n2'], ['h', 'co'], ['h', 'c2h4'], ['ch3', 'c2h4'],
#      ['hnc', 'hcn'], ['h', 'hcl'], ['oh', 'H2'], ['ch3', 'H2'], ['oh', 'CH4'], ['oh', 'NH3'], ['hcl', 'ch3'],
#      ['oh', 'C2H6'], ['f', 'H2'], ['O', 'CH4'], ['h', 'PH3'], ['h', 'oh'], ['h', 'H2S'], ['O', 'hcl'], ['NH2', 'ch3'],
#      ['NH2', 'C2H5'], ['C2H6', 'NH2'], ['NH2', 'CH4']]
# ref_l = [-64.91, -26.32, -103.28, -52.64, -32.16, -26.12, -20.32, -36.71, 3.69, -19.55, -40.02, -26.56, -15.06, -1.9,
#          -16.39, -3.11, -13.28, -10.32, -5.01, -16.84, -32.22, 5.44, -21.74, -2.32, -13.26, 0.42, -13.12, -9.56, -6.52,
#          -2.96]
wd_1 = '/home/guests/ychen/autorec/GMTKN55/FH51_benchmark'
#
# r = 30
# reactants = r_l[r-1]

wd_2 = f'/home/guests/ychen/autorec/GMTKN55/FH51_results/3_'

# end = ref_l[r-1]

reactants = ['C4H9CHO', 'H2O2']
if not os.path.isdir(f"{wd_2}"):
    os.makedirs(f'{wd_2}')
os.chdir(f'{wd_2}')
if os.path.isfile(f'{wd_2}/rec'):
    os.remove(f'{wd_2}/rec')

count = 1
sp_kcal_0 = 0

reatant_dirc = {}
reatant_dirc['1'] = reactants
for smi in reactants:
    # uhf = F.get_uhf(smi, wd_1)
    uhf =0
    F.dft_opt(smi=smi, m=uhf, wd_1=wd_1, wd_2=wd_2, c = 0)
    sp_kcal_0 = sp_kcal_0 + F.get_freeenergy(smi,wd_1,wd_2)
trige = 0
while count < 6:
    with open(f'{wd_2}/rec', 'a+') as f:
        f.write(f'count = {count}\n')

    if trige == 1:
        break
    for i in list(reatant_dirc.keys()):
        if len(i) == count:
            reactants_temp = reatant_dirc[i]

            intra_gap_list = []
            intra_rec_list = []
            for smi in reactants_temp:
                intra_gap_list.append(F.spr_intra(smi,wd_1,wd_2)[0])
                intra_rec_list.append(F.spr_intra(smi,wd_1,wd_2)[1])
            dis_dict = {}
            smi_count = 1
            for item in intra_gap_list:
                temp_dict = {item: intra_rec_list[intra_gap_list.index(item)]}
                dis_dict[str(smi_count)] = temp_dict
                smi_count = smi_count + 1

            dis_gap = list(dis_dict[list(dis_dict.keys())[0]].keys())[0]
            for a in list(dis_dict.keys()):
                b = list(dis_dict[a].keys())[1:]
                gap = list(dis_dict[a].keys())[0]
                if gap > dis_gap:
                    del dis_dict[a]
                else:
                    dis_gap = gap
            try:
                push_dict_1, push_dict_2 = F.spr_inter(reactants_temp,wd_1,wd_2)
                l1 = list(push_dict_1.keys())
                gap_1 = float('inf')
                gap_2 = float('inf')
                for item in l1:
                    a, b = item.split(",")
                    gap = list(push_dict_1[item].keys())[0]
                    if int(a) < int(b):
                        if gap < gap_1:
                            gap_1 = gap
                            push_smi_a1, push_smi_b1 = a, b
                    else:
                        if gap < gap_2:
                            gap_2 = gap
                            push_smi_a2, push_smi_b2 = a, b
                push_rec_1 = push_dict_1[f'{push_smi_a1},{push_smi_b1}'][gap_1]
                push_rec_2 = push_dict_1[f'{push_smi_a2},{push_smi_b2}'][gap_2]
                l2 = list(push_dict_2.keys())
                gap_3 = float('inf')
                gap_4 = float('inf')
                for item in l2:
                    a, b = item.split(",")
                    gap = list(push_dict_2[item].keys())[0]
                    if int(a) < int(b):
                        if gap < gap_3:
                            gap_3 = gap
                            push_smi_a3, push_smi_b3 = a, b
                    else:
                        if gap < gap_4:
                            gap_4 = gap
                            push_smi_a4, push_smi_b4 = a, b
                push_rec_3 = push_dict_2[f'{push_smi_a3},{push_smi_b3}'][gap_3]
                push_rec_4 = push_dict_2[f'{push_smi_a4},{push_smi_b4}'][gap_4]
            except:
                pass
            l_sp = {}
            reactants_1 = reactants_temp.copy()
            try:
                reactants_1, sp_kcal_1 = F.push(reactants_temp=reactants_1, push_rec=push_rec_1,
                                                push_smi_a=push_smi_a1,
                                                push_smi_b=push_smi_b1, count=count, num=str(i) + '1', wd_1=wd_1,
                                                wd_2=wd_2)
                l_sp['1'] = sp_kcal_1 - sp_kcal_0
                reatant_dirc[str(i) + '1'] = reactants_1
            except Exception as e:
                traceback.print_exc()
                pass

            reactants_2 = reactants_temp.copy()
            try:
                if push_smi_a2 == push_smi_b1 and push_smi_a1 == push_smi_b2:
                    if push_rec_1[0] == push_rec_2[0][::-1]:
                        reactants_2 = []
                reactants_2, sp_kcal_2 = F.push(reactants_temp=reactants_2, push_rec=push_rec_2,
                                                push_smi_a=push_smi_a2,
                                                push_smi_b=push_smi_b2, count=count, num=str(i) + '2', wd_1=wd_1,
                                                wd_2=wd_2)

                l_sp['2'] = sp_kcal_2 - sp_kcal_0
                reatant_dirc[str(i) + str(2)] = reactants_2
            except Exception as e:
                traceback.print_exc()
                pass


            reactants_3 = reactants_temp.copy()
            try:
                if push_smi_a3 == push_smi_a1 and push_smi_b3 == push_smi_b1:
                    if push_rec_3[0] == push_rec_1[0]:
                        reactants_3 = []
                if push_smi_a3 == push_smi_b2 and push_smi_b3 == push_smi_a2:
                    if push_rec_3[0] == push_rec_2[0][::-1]:
                        reactants_3 = []

                reactants_3, sp_kcal_3 = F.push(reactants_temp=reactants_3, push_rec=push_rec_3,
                                                push_smi_a=push_smi_a3,
                                                push_smi_b=push_smi_b3, count=count, num=str(i) + '3', wd_1=wd_1,
                                                wd_2=wd_2)
                l_sp['3'] = sp_kcal_3 - sp_kcal_0
                reatant_dirc[str(i) + str(3)] = reactants_3
            except Exception as e:
                traceback.print_exc()
                pass

            reactants_4 = reactants_temp.copy()
            try:
                if push_smi_a4 == push_smi_b1 and push_smi_b4 == push_smi_a1:
                    if push_rec_4[0] == push_rec_1[0][::-1]:
                        reactants_4 = []
                if push_smi_a4 == push_smi_a2 and push_smi_b4 == push_smi_b2:
                    if push_rec_4[0] == push_rec_2[0]:
                        reactants_4 = []
                if push_smi_a4 == push_smi_b3 and push_smi_b4 == push_smi_a3:
                    if push_rec_4[0] == push_rec_3[0][::-1]:
                        reactants_4 = []

                reactants_4, sp_kcal_4 = F.push(reactants_temp=reactants_4, push_rec=push_rec_4,
                                                push_smi_a=push_smi_a4,
                                                push_smi_b=push_smi_b4, count=count, num=str(i) + '4', wd_1=wd_1,
                                                wd_2=wd_2)
                l_sp['4'] = sp_kcal_4 - sp_kcal_0
                reatant_dirc[str(i) + str(4)] = reactants_2
            except Exception as e:
                print(e)
                pass

            # disso
            reactants_5 = reactants_temp.copy()
            try:
                reactants_5, sp_kcal_5 = F.dis(reactants_temp=reactants_5, dis_dict=dis_dict, count=count,
                                               num=str(i) + '5', wd_1=wd_1, wd_2=wd_2)
                l_sp['5'] = sp_kcal_5 - sp_kcal_0
                reatant_dirc[str(i) + str(5)] = reactants_5
            except Exception as e:
                traceback.print_exc()
                pass
            print(l_sp)
            with open(f'{wd_2}/rec', 'a+') as f:
                f.write(f'{l_sp}\n')
            if not l_sp:
                break


            for j in list(l_sp.keys()):
                if l_sp[j] > 500:
                    del l_sp[j]
                    del reatant_dirc[str(i) + str(j)]
                    if os.path.isdir(f'{wd_2}/{str(i) + str(j)}_push_{count}'):
                        shutil.rmtree(f'{wd_2}/{str(i) + str(j)}_push_{count}')
                    else:
                        shutil.rmtree(f'{wd_2}/{str(i) + str(j)}_disso_{count}')
            for j in list(l_sp.keys()):
                if l_sp[j] < -1000:
                    trige = 1
                    break

    with open(f'{wd_2}/rec', 'a+') as f:
        f.write(f'{reatant_dirc}\n')
    # print(reatant_dirc)
    # print(count)
    count = count + 1


