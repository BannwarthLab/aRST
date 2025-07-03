# updating reactants for this cycle
import shutil
from os.path import join,exists
import aRST.script.initialization as storage
from aRST.script.toolbox import get_random_folder_name,create_folder,get_strucinfo_from_path,write_xyz
from aRST.script.callsys import CallCrest


def updating_reacntans(reactants,job="depotonation"):
    # borrow gfn2 opt option for forcing gfn opt in refine_struc_info
    oldsetting = storage.setting.get_value("read.general.gfn2opt", False)
    storage.setting.set_value("read.general.gfn2opt", True)

    for struc in reactants:
        c = storage.allstruc.get_value(f"{struc}.struc_info.general.charge")
        constrain_l, ele = storage.allstruc.get_value(f"{struc}.struc_info.general.atomic_connection")
        if job=="depotonation":
            if not "H" in ele:
                continue

        E0 = storage.allstruc.get_value(f"{struc}.struc_info.energy.sp_dft")
        coord = storage.allstruc.get_value(f"{struc}.struc_info.general.coord")
        m = int(storage.allstruc.get_value(f"{struc}.struc_info.general.multiplicity"))

        tmpfolder = create_folder(join(storage.setting.bufferwd, get_random_folder_name(),job))
        # call depotonation
        write_xyz(coord, "coord")
        crest_parse = CallCrest(tmpfolder,
                         custom_settings={"charge": c, "uhf": m-1})
        if job=="depotonation":
            status = crest_parse.depotonation()
        elif job=="potonation":
            status = crest_parse.potonation()
        else:
            continue

        if status==0:
            updated_coord = crest_parse.read_update_coord(job=job)

            if updated_coord is not None:
                # reset AFO features
                storage.allstruc.set_value(f"{struc}.struc_info.reactivity.afo",None)
                storage.allstruc.set_value(f"{struc}.struc_info.reactivity.bo", None)
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

                storage.allstruc.refine_struc_info(updated_coord, c, m, storage.setting,struc,wd0=tmpfolder)
                struc_info = get_strucinfo_from_path(tmpfolder)
                storage.allstruc.set_value_from_struc_info(struc, struc_info)
                tarwd = join(storage.setting.strucwd, struc,job)
                if exists(tarwd):
                    shutil.rmtree(tarwd)
                shutil.move(tmpfolder,tarwd)
            else:
                print(f"The non-deprotonated/non-protonated form of the reactant {struc} is more energetically favorable->\n"
                      f"  keep the original reactant structure.")

    # reset gfn2 opt option
    storage.setting.set_value("read.general.gfn2opt", oldsetting)