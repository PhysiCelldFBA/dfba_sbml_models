#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
import cobra
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from cobra.flux_analysis import find_blocked_reactions

tolerance = 1e-9


# In[2]:


def update_exchange_dict(model, summary, exchanges_dict, tolerance=1e-7):
    print(f"Secretion fluxes for {model.id}")
    for r in summary.secretion_flux.index:
        flux = -1 * summary.secretion_flux.loc[r, 'flux']
        ub = model.reactions.get_by_id(r).upper_bound
        if abs(flux) < tolerance:
            continue
        print(f" - {r} {flux} ({ub})")
        if r not in exchanges_dict:
            exchanges_dict[r] = [0.0, 1000.0]
        exchanges_dict[r][1] = max(exchanges_dict[r][1], ub)
    
    print()
    print(f"Uptake fluxes for {model.id}")
    for r in summary.uptake_flux.index:
        flux = -1 * summary.uptake_flux.loc[r, 'flux']
        if abs(flux) < tolerance:
            continue
       
        lb = model.reactions.get_by_id(r).lower_bound
        print(f" - {r} {flux} ({lb})")
        if r not in exchanges_dict:
            exchanges_dict[r] = [0.0, 1000.0]
        exchanges_dict[r][0] = min(exchanges_dict[r][0], lb)

    return exchanges_dict


# In[3]:


# Fixing IDs for reaction changinx ex prefix to EX
cb_model = read_sbml_model("models/iCB925.xml")
cb_model.tolerance = tolerance
cb_model.id = "iCB925"
for r in cb_model.boundary:
    r.id = "EX" + r.id[2:]
    r.name = "EX" + r.name[2:]

cb_model.repair()
for r in cb_model.boundary:
    if r.id.startswith("EX"):
        continue
    print(r.id)

print(f"Exchange reaction renamed correctly for {cb_model.id}")
cb_sol = cb_model.optimize()
cb_summary = cb_model.summary()
cb_summary


# In[4]:


mb_model = read_sbml_model("models/iAF692.xml")
mb_model.tolerance = tolerance
for r in sorted(mb_model.boundary, key=lambda x: x.id):
    r_id = re.sub("__", "_", r.id)
    r.id = r_id

mb_model.repair()
for r in mb_model.boundary:
    if r.id.startswith("EX"):
        continue
    print(r.id)

print(f"Exchange reaction renamed correctly for {mb_model.id}")
print()
print("Fixing bounds to remove meoh from the media and cysteine and enable H2")
mb_model.reactions.get_by_id("EX_meoh_e").lower_bound = 0
mb_model.reactions.get_by_id("EX_h2_e").lower_bound = -20
mb_model.reactions.get_by_id("EX_cys_L_e").lower_bound = 0
print()

mb_sol = mb_model.optimize()
mb_summary = mb_model.summary()
mb_summary


# In[5]:


exchanges_dict = {}
exchanges_dict = update_exchange_dict(cb_model, cb_summary, exchanges_dict)
exchanges_dict = update_exchange_dict(mb_model, mb_summary, exchanges_dict)

print("\n======================================================\n")

for r in sorted(exchanges_dict):
    global_bounds = tuple(exchanges_dict[r])
    print(f" - {r}")
    
    if r in cb_model.reactions:
        rxn = cb_model.reactions.get_by_id(r)
        print(f"\tCb bounds", rxn.bounds)
    else:
        print(f"\tCb bounds NONE")
    if r in mb_model.reactions:
        rxn = mb_model.reactions.get_by_id(r)
        print(f"\tMb bounds", rxn.bounds)
    else:
        print(f"\tMb bounds NONE")
        
    print(f"\tGlobal bounds", global_bounds)


# In[6]:


print(f"Updating {cb_model.id} model")
new_cb_model = cb_model.copy()
for r in new_cb_model.boundary:
    if r.id not in exchanges_dict:
        r.bounds = (0.0, 1000.0)
        continue
    else:
        lb, ub = exchanges_dict[r.id]
        r.lower_bound = lb
        r.upper_bound = ub

cb_sol = cb_model.optimize()
new_cb_sol = new_cb_model.optimize()
print(new_cb_sol.objective_value, cb_sol.objective_value)
ratio = abs(new_cb_sol.objective_value - cb_sol.objective_value) / cb_sol.objective_value * 100
if ratio > 10:
    print("WARNING")


cb_blocked = find_blocked_reactions(new_cb_model, processes=12, zero_cutoff=tolerance)
cb_gap_metabolties = []
for m in new_cb_model.metabolites:
    m_rxns = {r.id for r in m.reactions}
    if len( m_rxns - set(cb_blocked) ) == 0:
        cb_gap_metabolties.append(m.id)

n_reaction = len(new_cb_model.reactions)
n_blocked = len(cb_blocked)
n_metabolites = len(new_cb_model.metabolites)
n_gap = len(cb_gap_metabolties)

print(f"{new_cb_model.id}: ")
print(f"\tReactions total/blocked: {n_reaction}/{n_blocked}")
print(f"\tMetabolites total/blocked: {n_metabolites}/{n_gap}")
print()
print("Trimming model by removing blocked reactions and gap metabolites")
print()

new_cb_model.remove_reactions(cb_blocked)
to_remove = [new_cb_model.metabolites.get_by_id(m) for m in cb_gap_metabolties]
new_cb_model.remove_metabolites(to_remove)
new_cb_sol = new_cb_model.optimize()

print(new_cb_sol.objective_value, cb_sol.objective_value)

ratio = abs(new_cb_sol.objective_value - cb_sol.objective_value) / cb_sol.objective_value * 100
if ratio > 10:
    print("WARNING")
    new_cb_model.summary()
else:
    fname = f"models/{new_cb_model.id}_trimmed.xml"
    print(f"Saving {new_cb_model.id} trimmed into {fname}")
    write_sbml_model(new_cb_model, fname)
    n_reaction = len(new_cb_model.reactions)
    n_metabolites = len(new_cb_model.metabolites)
    print(f"{new_cb_model.id}: ")
    print(f"\tReactions total: {n_reaction}")
    print(f"\tMetabolites total: {n_metabolites}")


# In[7]:


print(f"Updating {mb_model.id} model")
new_mb_model = mb_model.copy()
for r in new_mb_model.boundary:
    if r.id not in exchanges_dict:
        r.bounds = (0.0, 1000.0)
        continue
    else:
        lb, ub = exchanges_dict[r.id]
        r.lower_bound = lb
        r.upper_bound = ub

new_mb_sol = new_mb_model.optimize()
print(new_mb_sol.objective_value, mb_sol.objective_value)
ratio = abs(new_mb_sol.objective_value - mb_sol.objective_value) / mb_sol.objective_value * 100
if ratio > 10:
    print("WARNING")


mb_blocked = find_blocked_reactions(new_mb_model, processes=12, zero_cutoff=1e-9)
mb_gap_metabolties = []
for m in new_mb_model.metabolites:
    m_rxns = {r.id for r in m.reactions}
    if len( m_rxns - set(mb_blocked) ) == 0:
        mb_gap_metabolties.append(m.id)

n_reaction = len(new_mb_model.reactions)
n_blocked = len(mb_blocked)
n_metabolites = len(new_mb_model.metabolites)
n_gap = len(mb_gap_metabolties)
print(f"{new_mb_model.id}: ")
print(f"\tReactions total/blocked: {n_reaction}/{n_blocked}")
print(f"\tMetabolites total/blocked: {n_metabolites}/{n_gap}")
print()
print("Trimming model by removing blocked reactions and gap metabolites")
print()
new_mb_model.remove_reactions(mb_blocked)
to_remove = [new_mb_model.metabolites.get_by_id(m) for m in mb_gap_metabolties]
new_mb_model.remove_metabolites(to_remove)
new_mb_sol = new_mb_model.optimize()

print(new_mb_sol.objective_value, mb_sol.objective_value)

ratio = abs(new_mb_sol.objective_value - mb_sol.objective_value) / mb_sol.objective_value * 100
if ratio > 10:
    print("WARNING")
    new_mb_model.summary()
else:
    n_reaction = len(new_mb_model.reactions)
    n_metabolites = len(new_mb_model.metabolites)
    print(f"{new_mb_model.id}: ")
    print(f"\tReactions total: {n_reaction}")
    print(f"\tMetabolites total: {n_metabolites}")
    new_mb_sol = new_mb_model.optimize()
    print(f"\tOptimal solution: {new_mb_sol.objective_value}")
    fname = f"models/{new_mb_model.id}_trimmed.xml"
    print(f"Saving {new_mb_model.id} trimmed into {fname}")
    write_sbml_model(new_mb_model, fname)


# In[8]:


for model_id in ["iCB925", "iAF692"]: 
    fname = f"models/{model_id}_trimmed.xml"
    print(f"Reading {fname}")
    model = read_sbml_model(fname)
    n_reaction = len(model.reactions)
    n_metabolites = len(model.metabolites)
    print(f"{model.id}: ")
    print(f"\tReactions total: {n_reaction}")
    print(f"\tMetabolites total: {n_metabolites}")
    sol = model.optimize()
    print(f"\tOptimal solution: {sol.objective_value}")

