# updating reactants for this cycle
import shutil
from os.path import join,exists
import aRST.script.initialization as storage
from aRST.script.toolbox import get_random_folder_name,create_folder,get_strucinfo_from_path,write_xyz
from aRST.script.callsys import CallCrest


def updating_reacntans(reactants,job="depotonation"):
    # borrow gfn2 opt option for forcing gfn opt in refine_struc_info
    oldsetting = storage.setting.read.general.gfn2opt
    storage.setting.read.general.gfn2opt=True

    for struc in reactants:
        s = storage.allstruc.get(struc)
        c = s.struc_info.general.charge
        constrain_l = s.struc_info.general.atomic_connection
        ele_l = s.struc_info.general.ele_l
        if job=="depotonation":
            if not "H" in ele_l:
                continue

        s = storage.allstruc.get(struc)
        E0 = s.struc_info.energy.sp
        coordinate = s.struc_info.general.coordinate
        ele_l = s.struc_info.general.ele_l
        m = s.struc_info.general.multiplicity

        tmpfolder = create_folder(join(storage.setting.searching_setting.general.bufferwd, get_random_folder_name(),job))
        # call depotonation
        write_xyz(coordinate,ele_l, "coord")
        crest_parse = CallCrest(tmpfolder,
                         custom_settings={"charge": c, "uhf": m-1})
        if job=="depotonation":
            status = crest_parse.depotonation()
        elif job=="potonation":
            status = crest_parse.potonation()
        else:
            continue

        if status==0:
            updated_coordinate,updated_ele_l = crest_parse.read_update_coord(job=job)

            if updated_coordinate is not None:
                # reset AFO features
                s = storage.allstruc.get(struc)
                s.struc_info.reactivity.afo = None
                s.struc_info.reactivity.bo = None
                if job == "depotonation":
                    print(
                        f"The deprotonated form of the reactant {struc} is more energetically favorable->\n"
                        f"  updated to reactant")
                    c -= 1
                else:
                    print(
                        f"The protonated form of the reactant {struc} is more energetically favorable->\n"
                        f"  updated to reactant")
                    c += 1

                storage.allstruc.refine_struc_info(updated_coordinate,updated_ele_l, c, m, storage.setting,strucid=struc,wd0=tmpfolder)
                struc_info = get_strucinfo_from_path(tmpfolder)
                storage.allstruc.set_value_from_struc_info(struc, struc_info)
                tarwd = join(storage.setting.searching_setting.general.strucwd, struc,job)
                if exists(tarwd):
                    shutil.rmtree(tarwd)
                shutil.move(tmpfolder,tarwd)
            else:
                print(f"The non-deprotonated/non-protonated form of the reactant {struc} is more energetically favorable->\n"
                      f"  keep the original reactant structure.")

    # reset gfn2 opt option
    storage.setting.read.general.gfn2opt=oldsetting