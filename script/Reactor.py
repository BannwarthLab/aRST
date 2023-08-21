# def do_FULLreaction(wd, react_l, rec, **kwargs):
#     xtbwd = kwargs.get('xtbwd')
#     alpb = kwargs.get('alpb', False)
#     level = kwargs.get('dft_level', 'r2SCAN-3c')
#     energy_lim = kwargs.get('elim', 200)
#     # ml features
#     if alpb:
#         ml(rec, react_l, wd, xtbwd, alpb)
#     else:
#         ml(rec, react_l, wd, xtbwd)
#     rec = collect_ml(rec, react_l, wd)
#
#     # Asoociation
#     if alpb:
#         reactants_inter_dict = push(react_l, rec, wd, alpb, kwargs=kwargs)
#     else:
#         reactants_inter_dict = push(react_l, rec, wd, kwargs=kwargs)
#     ninter = len(reactants_inter_dict.keys())
#
#     # Dissociation
#     if alpb:
#         reactants_dis_dict = pull(react_l, rec, wd, alpb, kw={'association_count': ninter})
#     else:
#         reactants_dis_dict = pull(react_l, rec, wd, kw={'association_count': ninter})
#
#     # Subprocessing after scan
#     if alpb:
#         rec = do_subprocess(wd, react_l, reactants_inter_dict, reactants_dis_dict, rec, level,energy_lim, alpb)
#
#     else:
#         rec = do_subprocess(wd, react_l, reactants_inter_dict, reactants_dis_dict, rec, level,energy_lim)
#
#     return rec