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


def model_trimming(model, eps_tolerance=1e-9, bound_scaling = 1.0, cap_internal_bounds=False):
    """Remove the reactions that are blocked in the medium of `model`, the metabolites left without
    reactions and the unconnected genes.

    The bounds of the internal reactions are the ones of `model`. With cap_internal_bounds=True
    (the behaviour of the first release of the trimmed models) every internal reaction is instead
    bounded to its flux range in this medium, multiplied by bound_scaling; these caps restrict the
    model in any other medium (iMG746: optimum 15.7% lower, 314 upper bounds cut from 1000 to below 1)
    and are not used. The exchange reactions always take their flux range in this medium.
    """
    cmodel = model.copy()

    fva_df = flux_variability_analysis(cmodel, processes=10, fraction_of_optimum=0.0)

    blocked_reaction = fva_df.index[ fva_df.abs().max(axis=1) < eps_tolerance ]
    
    mask = fva_df['minimum'].abs() < eps_tolerance
    fva_df.loc[mask, 'minimum'] = 0.0
    
    mask = fva_df['maximum'].abs() < eps_tolerance
    fva_df.loc[mask, 'maximum'] = 0.0
    
    mask = fva_df['minimum'] > eps_tolerance
    fva_df.loc[mask, 'minimum'] = 0.0
    
    maks = fva_df['maximum'] > -eps_tolerance
    fva_df.loc[mask, 'maximum'] = 0.0

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
    for rxn in cmodel.reactions:
        is_exchange = rxn.id.startswith("EX") or rxn.id.startswith("ex")
        if is_exchange or cap_internal_bounds:
            rxn.lower_bound = float(fva_df.loc[rxn.id, 'minimum'])
            rxn.upper_bound = float(fva_df.loc[rxn.id, 'maximum'])

    if cap_internal_bounds:
        for r in cmodel.reactions:
            if r.id.startswith("EX") or r.id.startswith("ex"):
                continue
            if r.upper_bound > eps_tolerance:
                r.upper_bound *= bound_scaling
            if r.lower_bound < -eps_tolerance:
                r.lower_bound *= bound_scaling

    return cmodel

    