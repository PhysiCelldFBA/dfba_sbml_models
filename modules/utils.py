import re
import numpy as np
from cobra.manipulation.delete import remove_genes
from cobra.flux_analysis import flux_variability_analysis


def correct_ids(model, pattern, fix):
    total = 0
    for r in model.reactions:
        if re.search(pattern, r.id):
           r.id = re.sub(pattern, fix, r.id)
           total += 1
    print(f"N of corrected reaction ids {total}")
    model.repair()
    for m in model.metabolites:
        if re.search(pattern, m.id):
           m.id = re.sub(pattern, fix, m.id)
           total += 1

    print(f"N of corrected metabolites ids {total}")
    model.repair()
    return model


def trimm_model(model, eps_tolerance=1e-8, bound_scaling = 1.2):
    cmodel = model.copy()

    fva_dict = flux_variability_analysis(cmodel, processes=10, fraction_of_optimum=0.0)
    
    blocked_reaction = []
    new_bounds = {}
    
    for r in fva_dict.index:
        lb = fva_dict.loc[r, 'minimum']
        ub = fva_dict.loc[r, 'maximum']
        if max(np.abs([lb, ub])) < eps_tolerance:
            blocked_reaction.append(r)
        else:
            if abs(ub) < eps_tolerance:
                ub = 0.0
            if abs(lb) < eps_tolerance:
                lb = 0.0
            new_bounds[r] = (lb, ub)

    print(f"Removing blocked {len(blocked_reaction)} reactions")
    cmodel.remove_reactions(blocked_reaction)
    
    gap_metabolites = []
    for m in cmodel.metabolites:
        R = {r.id for r in m.reactions}
        if len(R - set(blocked_reaction)) == 0:
            gap_metabolites.append(m)
    
    print(f"Removing {len(gap_metabolites)} gap metabolites")
    cmodel.remove_metabolites(gap_metabolites)
    
    unconnected_genes = []
    for g in cmodel.genes:
        if len(g.reactions) == 0:
            unconnected_genes.append(g)
    
    print(f"Removing {len(unconnected_genes)} unconnected genes")
    remove_genes(cmodel, unconnected_genes)
    
    print("Setting bounds")
    for r,(lb,ub) in new_bounds.items():
        rxn = cmodel.reactions.get_by_id(r)
        rxn.upper_bound = float(ub)
        rxn.lower_bound = float(lb)
    
    
    for r in cmodel.reactions:
        if r.id.startswith("EX") or r.id.startswith("ex"):
            continue
        if r.upper_bound > eps_tolerance:
            r.upper_bound *= bound_scaling 
        if r.lower_bound < -eps_tolerance:
            r.lower_bound *= bound_scaling

    return cmodel

    