import os
import shutil
import billiard as multiprocessing
from Assign_reaction import do_assign_reactor, do_reaction
from operation import write_rec
from read import update_dict
from Subprocess import elvauate_reaction


def reactionnetwork(rec, wd, **kwargs):
    # rec: Record all reactants info: name, coord, at_features, sp
    # rec_valid: Record non-repeat reactants info: name, coord, at_features, sp
    # reaction_rec: Record valid reactions paths {path_name: [[R],[P],[delt_E]}
    # wd : Working folder
    # cycle: max simulation cycles
    if kwargs.get('kwargs'):
        kwargs = kwargs.get('kwargs')
    cycle = kwargs.get('cycle', 1)
    reactor = kwargs.get('reactor_like', False)
    energy_lim = kwargs.get('elim', 200)

    reaction_rec = {}
    rec_valid = rec
    step = 0
    while step < cycle:
        if reactor:
            react_l = [k for k in list(rec_valid.keys())]  # all species appeared in system
            print(f"nr of current species: {len(react_l)}")
            print(f"rec_l for step {step}: {react_l} ")
            reactants_tmp = do_assign_reactor(react_l, rec, reaction_rec)  # random select two molecules
            print("reactants = ", reactants_tmp)

            reactants_total_dict, rec_tmp = do_reaction(wd, rec,
                                                        reactants_tmp, step + 1, kwargs=kwargs)  # 1 asso + 1 disso
            rec = update_dict(rec, rec_tmp)

        else:
            path_l = list(set([k.split('/')[0] for k in list(rec.keys())
                               if len(k.split('/')[0]) == step + 1
                               and os.path.isdir(f'{wd}/{k.split("/")[0]}')]))  # reaction path for this cycle

            # individually sim each path
            with multiprocessing.Pool(len(path_l)) as pool:
                results = []
                for i in range(len(path_l)):
                    path_name = path_l[i]
                    reactants_tmp = [k for k in list(rec.keys())
                                     if k.split("/")[0] == path_name]  # record reactants for this cycle
                    print(f"path: {path_name}; reactants: {reactants_tmp}")
                    results.append(
                        pool.apply_async(do_reaction, (wd, rec, reactants_tmp), kwargs))
                pool.close()
                pool.join()
            reactants_total_dict = {}
            for res in results:
                reactants_total_dict = update_dict(reactants_total_dict, res.get()[0])
                rec = update_dict(rec, res.get()[1])

        if reactants_total_dict:
            # collect energy and evaluate reaction
            path_l = list(reactants_total_dict.keys())
            del_l = []
            for path_i in path_l:
                reaction_rec_tmp, del_l_tmp = elvauate_reaction(path_i, wd, rec, energy_lim, reactants_total_dict,
                                                                reaction_rec, reactor)
                reaction_rec = update_dict(reaction_rec, reaction_rec_tmp)
                del_l = del_l + del_l_tmp

            print("path_l: ", path_l)
            print("del_l: ", list(set(del_l)))

            # kill paths that disqualified
            for i in list(set(del_l)):
                if os.path.isdir(f'{wd}/{i}'):
                    shutil.rmtree(f'{wd}/{i}', ignore_errors=True)
                for j in list(rec.keys()):
                    if j.split('/')[0] == i:
                        del rec[j]
                del reaction_rec[i]
                path_l.remove(i)

            # update rec_valid
            rec_valid = {}
            tmp = []
            for i in list(rec.keys()):
                if 'mb' in rec[i].keys():
                    mb = rec[i]['mb']
                    if mb not in tmp:
                        tmp.append(mb)
                        rec_valid[i] = {}
                        rec_valid[i] = rec[i]
            # write document
            if reactor:
                if path_l:  # step +1 only if there is association [s for s in path_l if "_1" in s]
                    write_rec(rec_valid, wd, 'molecules')
                    write_rec(reaction_rec, wd, 'reaction')
                    step += 1
                else:
                    continue
            else:
                if path_l:
                    write_rec(rec_valid, wd, 'molecules')
                    write_rec(reaction_rec, wd, 'reaction')
                    step += 1
                else:
                    break
        else:
            break
    else:
        return rec
