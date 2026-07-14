#!/usr/bin/env python3
"""Orchestrate the repository model preparation workflow.

This script replaces the old monolithic trimming script with a simple
entry-point that prepares each curated model in turn.
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Tuple

import cobra
import numpy as np
from cobra.io import read_sbml_model, validate_sbml_model, write_sbml_model
from cobra.flux_analysis import flux_variability_analysis
from cobra.manipulation.delete import remove_genes

try:
    from scipy.io import loadmat
    from scipy.sparse import csc_matrix
except ImportError:  # pragma: no cover - optional dependency for MCF7 conversion
    loadmat = None
    csc_matrix = None

REPO_ROOT = Path(__file__).resolve().parent
MODELS_DIR = REPO_ROOT / "models"
MODULES_DIR = REPO_ROOT / "modules"

if str(MODULES_DIR) not in sys.path:
    sys.path.insert(0, str(MODULES_DIR))

from utils import correct_ids, model_trimming


def prepare_clostridium() -> None:
    """Prepare the Clostridium beijerinckii model."""
    print("Preparing Clostridium beijerinckii model")

    original_model_path = MODELS_DIR / "iCB925.xml"
    final_model_path = MODELS_DIR / "C_beijerinckii_model_iCB925_trimmed.xml"

    cb_model = read_sbml_model(str(original_model_path))
    cb_model = correct_ids(cb_model, "ex_", "EX_")

    cb_model.reactions.EX_glc_D_e.lower_bound = -7.5
    cb_model.reactions.ATPM.lower_bound = 0

    solution = cb_model.optimize()
    print(f"Initial biomass objective: {solution.objective_value}")

    trimmed_model = model_trimming(cb_model, eps_tolerance=1e-9, bound_scaling=10.0)
    for reaction in trimmed_model.exchanges:
        if reaction.upper_bound > 0:
            reaction.upper_bound = 100

    trimmed_model.reactions.biomass.upper_bound = 0.35
    for reaction_id in ["EX_h2o_e", "EX_h_e", "EX_nh4_e", "EX_pi_e", "EX_so4_e"]:
        trimmed_model.reactions.get_by_id(reaction_id).lower_bound = -100

    solution_trimmed = trimmed_model.optimize()
    print(f"Trimmed biomass objective: {solution_trimmed.objective_value}")
    assert abs(solution.objective_value - solution_trimmed.objective_value) < 1e-9

    trimmed_model.notes = {}
    for metabolite in trimmed_model.metabolites:
        metabolite.notes = {}

    write_sbml_model(trimmed_model, str(final_model_path))
    print(f"Wrote {final_model_path}")


def prepare_methanosarcina() -> None:
    """Prepare the Methanosarcina barkeri model."""
    print("Preparing Methanosarcina barkeri model")

    original_model_path = MODELS_DIR / "M_barkeri_model_iMG746" / "mba_model.xml"
    intermediate_model_path = MODELS_DIR / "M_barkeri_model_iMG746.xml"
    final_model_path = MODELS_DIR / "M_barkeri_model_iMG746_trimmed.xml"

    mb_model = read_sbml_model(str(intermediate_model_path))
    mb_model, errors = validate_sbml_model(str(original_model_path))

    assert len(errors["SBML_FATAL"]) == 0
    assert len(errors["SBML_ERROR"]) == 0
    assert len(errors["SBML_SCHEMA_ERROR"]) == 0

    replace_dict = {
        "_LPAREN_": "_",
        "_RPAREN_": "",
        "_DASH_": "_",
        "_unknown_": "_",
        "unknown_": "",
    }

    for pattern, fix in replace_dict.items():
        mb_model = correct_ids(mb_model, pattern, fix)

    write_sbml_model(mb_model, str(intermediate_model_path))
    mb_model = read_sbml_model(str(intermediate_model_path))

    mb_model.reactions.ATPM.lower_bound = 0
    mb_model.reactions.EX_h2_e.lower_bound = -20

    mb_model.reactions.EX_4abz_e.lower_bound = -1
    mb_model.reactions.EX_co2_e.lower_bound = -100
    mb_model.reactions.EX_cobalt2_e.lower_bound = -1
    mb_model.reactions.EX_h2s_e.lower_bound = -10
    mb_model.reactions.EX_nac_e.lower_bound = -10
    mb_model.reactions.EX_nh4_e.lower_bound = -10
    mb_model.reactions.EX_ni2_e.lower_bound = -10
    mb_model.reactions.EX_pi_e.lower_bound = -10
    mb_model.reactions.EX_meoh_e.lower_bound = 0
    mb_model.reactions.EX_cys_L_e.lower_bound = 0

    solution_co2 = mb_model.optimize()
    print(f"Growth rate with CO2 as carbon source: {solution_co2.objective_value:.5f}")

    mb_model.reactions.get_by_id("EX_ac_e").lower_bound = -3
    solution_acetate = mb_model.optimize()
    print(f"Growth rate with acetate as carbon source: {solution_acetate.objective_value:.5f}")

    trimmed_model = model_trimming(mb_model, bound_scaling=10)
    for reaction in trimmed_model.exchanges:
        if reaction.upper_bound > 0:
            reaction.upper_bound = 100

    trimmed_model.reactions.Mb_biomass_65.upper_bound = 0.05
    trimmed_model.reactions.EX_4abz_e.lower_bound = -1
    trimmed_model.reactions.EX_co2_e.lower_bound = -100
    trimmed_model.reactions.EX_cobalt2_e.lower_bound = -1
    trimmed_model.reactions.EX_h2s_e.lower_bound = -10
    trimmed_model.reactions.EX_nac_e.lower_bound = -10
    trimmed_model.reactions.EX_nh4_e.lower_bound = -10
    trimmed_model.reactions.EX_ni2_e.lower_bound = -10
    trimmed_model.reactions.EX_pi_e.lower_bound = -10

    trimmed_model.reactions.get_by_id("EX_ac_e").lower_bound = 0
    new_solution_co2 = trimmed_model.optimize()
    diff = solution_co2.objective_value - new_solution_co2.objective_value
    if abs(diff) > 1e-4:
        print(f"Warning: CO2 objective shifted by {diff:.6f}")

    trimmed_model.reactions.get_by_id("EX_ac_e").lower_bound = -3
    new_solution_ac = trimmed_model.optimize()
    diff = solution_acetate.objective_value - new_solution_ac.objective_value
    if abs(diff) > 1e-4:
        print(f"Warning: acetate objective shifted by {diff:.6f}")

    for metabolite in trimmed_model.metabolites:
        metabolite.notes = {}
    write_sbml_model(trimmed_model, str(final_model_path))
    print(f"Wrote {final_model_path}")


def trim_models() -> None:
    """Sanity-check and summarize the curated trimmed models."""
    print("Reviewing trimmed models")
    for path in [
        MODELS_DIR / "C_beijerinckii_model_iCB925_trimmed.xml",
        MODELS_DIR / "M_barkeri_model_iMG746_trimmed.xml",
        MODELS_DIR / "iML1515_trimmed_glc_fva99.xml",
        MODELS_DIR / "iML1515_trimmed_acetate_fva99.xml",
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Expected trimmed model missing: {path}")
        model = read_sbml_model(str(path))
        print(f"- {model.id}: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")


def prepare_iml1515() -> None:
    """Prepare iML1515 under simple aerobic glucose and acetate minimal media.

    This creates two trimmed artifacts that preserve near-optimal growth
    using FVA at 99%% of the optimum.
    """
    print("Preparing iML1515 model (glucose/acetate minimal media)")

    original_model_path = MODELS_DIR / "iML1515.xml"
    out_glc = MODELS_DIR / "iML1515_trimmed_glc_fva99.xml"
    out_ac = MODELS_DIR / "iML1515_trimmed_acetate_fva99.xml"

    if not original_model_path.exists():
        raise FileNotFoundError(f"Missing iML1515 model: {original_model_path}")

    model = read_sbml_model(str(original_model_path))

    # Ensure the canonical biomass reaction is selected as objective when present
    biomass_candidates = [r for r in model.reactions if r.id.startswith("BIOMASS")]
    if biomass_candidates:
        for rxn in model.reactions:
            rxn.objective_coefficient = 0
        biomass_candidates[0].objective_coefficient = 1
        print(f"Set objective to: {biomass_candidates[0].id}")

    common_ions = [
        "EX_nh4_e",
        "EX_pi_e",
        "EX_so4_e",
        "EX_k_e",
        "EX_na1_e",
        "EX_h2o_e",
        "EX_h_e",
        "EX_ca2_e",
        "EX_mg2_e",
    ]

    def _trim_preserve_optimal(m_in, fraction=0.99, eps=1e-9, bound_scaling=10.0):
        m = m_in.copy()
        fva = flux_variability_analysis(m, processes=4, fraction_of_optimum=fraction)
        blocked = fva.index[fva.abs().max(axis=1) < eps]
        if len(blocked) > 0:
            m.remove_reactions(blocked)
        gap_mets = [met for met in m.metabolites if all(r.id in blocked for r in met.reactions)]
        if gap_mets:
            m.remove_metabolites(gap_mets)
        unconnected = [g for g in m.genes if len(g.reactions) == 0]
        if unconnected:
            remove_genes(m, unconnected)
        # widen bounds to be conservative
        for rxn in m.reactions:
            if rxn.upper_bound > 0:
                rxn.upper_bound *= bound_scaling
            if rxn.lower_bound < 0:
                rxn.lower_bound *= bound_scaling
        return m

    # Glucose aerobic minimal
    mg = model.copy()
    if "EX_glc__D_e" in mg.reactions:
        mg.reactions.EX_glc__D_e.lower_bound = -10.0
    if "EX_o2_e" in mg.reactions:
        mg.reactions.EX_o2_e.lower_bound = -20.0
    for rid in common_ions:
        if rid in mg.reactions:
            mg.reactions.get_by_id(rid).lower_bound = -100.0

    print("Glucose pre-trim objective:", mg.optimize().objective_value)
    trimmed_g = _trim_preserve_optimal(mg, fraction=0.99, bound_scaling=10.0)
    for r in trimmed_g.exchanges:
        if r.upper_bound > 0:
            r.upper_bound = 100
    write_sbml_model(trimmed_g, str(out_glc))
    print(f"Wrote {out_glc}")
    print("Glucose trimmed objective:", trimmed_g.optimize().objective_value)

    # Acetate aerobic minimal
    ma = model.copy()
    if "EX_ac_e" in ma.reactions:
        ma.reactions.EX_ac_e.lower_bound = -10.0
    if "EX_o2_e" in ma.reactions:
        ma.reactions.EX_o2_e.lower_bound = -20.0
    for rid in common_ions:
        if rid in ma.reactions:
            ma.reactions.get_by_id(rid).lower_bound = -100.0

    print("Acetate pre-trim objective:", ma.optimize().objective_value)
    trimmed_a = _trim_preserve_optimal(ma, fraction=0.99, bound_scaling=10.0)
    for r in trimmed_a.exchanges:
        if r.upper_bound > 0:
            r.upper_bound = 100
    write_sbml_model(trimmed_a, str(out_ac))
    print(f"Wrote {out_ac}")
    print("Acetate trimmed objective:", trimmed_a.optimize().objective_value)


def _normalize_gene_id(gid: str) -> str:
    """Strip trailing .digits from gene IDs (e.g. '4830.1' -> '4830')."""
    if gid is None:
        return ""
    gid = str(gid).strip()
    return re.sub(r"\.(\d+)$", r"\1", gid)


def _normalize_gpr(gpr_str: str) -> str:
    """Clean GPR strings and remove numeric suffixes."""
    if gpr_str is None:
        return ""
    s = str(gpr_str).strip()
    if not s or s in {"[]", "nan", "None"}:
        return ""
    s = re.sub(r"(\d+)\.\d+", r"\1", s)
    return s.replace("{", "").replace("}", "")


def _transform_matlab_met_id(raw_met: str) -> Tuple[str, str]:
    """Transform 'met[c]' -> ('met_c', 'c')."""
    if raw_met is None:
        return "unknown", "c"
    raw = str(raw_met)
    match = re.match(r"^(.*)\[([A-Za-z0-9_]+)\]$", raw)
    if match:
        base, comp = match.group(1), match.group(2)
        return f"{base}_{comp}", comp
    return raw, "c"


def _cobramat2cobrapy(cobramat, name: str | None = None) -> cobra.Model:
    """Convert a MATLAB COBRA struct into a COBRApy model."""
    model_id = name or str(getattr(cobramat, "description", "unnamed_model"))
    model = cobra.Model(model_id)

    rxns = np.array(cobramat.rxns, dtype=str)
    mets = np.array(cobramat.mets, dtype=str)
    stoich = csc_matrix(cobramat.S)
    lb = np.array(cobramat.lb, dtype=float)
    ub = np.array(cobramat.ub, dtype=float)
    cvec = np.array(cobramat.c, dtype=float)
    b = np.array(getattr(cobramat, "b", np.zeros(len(mets))), dtype=float)

    met_formulas = getattr(cobramat, "metFormulas", None)
    met_names = getattr(cobramat, "metNames", None)
    gr_rules = getattr(cobramat, "grRules", None)
    genes = getattr(cobramat, "genes", None)
    rxn_names = getattr(cobramat, "rxnNames", None)
    sub_systems = getattr(cobramat, "subSystems", None)

    met_objs: Dict[str, cobra.Metabolite] = {}
    seen_ids = defaultdict(int)
    existing_met_ids = set()

    for i, raw_met in enumerate(mets):
        new_id, comp = _transform_matlab_met_id(raw_met)
        if new_id in existing_met_ids:
            seen_ids[new_id] += 1
            unique_id = f"{new_id}__dup{seen_ids[new_id]}"
        else:
            unique_id = new_id
            seen_ids[new_id] = 0
        existing_met_ids.add(unique_id)

        met = cobra.Metabolite(unique_id, compartment=comp)
        if met_formulas is not None and i < len(met_formulas):
            formula = str(met_formulas[i]).strip()
            if formula and formula not in {"[]", "nan", "None"}:
                met.formula = formula
        if met_names is not None and i < len(met_names):
            name_str = str(met_names[i]).strip()
            if name_str and name_str not in {"[]", "nan", "None"}:
                met.name = name_str
        met_objs[str(raw_met)] = met

    model.add_metabolites(list(met_objs.values()))

    coo = stoich.tocoo()
    col_to_entries = defaultdict(list)
    for row_idx, col_idx, value in zip(coo.row, coo.col, coo.data):
        col_to_entries[col_idx].append((row_idx, value))

    for j, rxn_id in enumerate(rxns):
        rxn = cobra.Reaction(rxn_id)
        if rxn_names is not None and j < len(rxn_names):
            rxn.name = str(rxn_names[j]).strip()
        if sub_systems is not None and j < len(sub_systems):
            rxn.subsystem = str(sub_systems[j]).strip()

        rxn.lower_bound = float(lb[j])
        rxn.upper_bound = float(ub[j])
        model.add_reactions([rxn])

        entries = col_to_entries.get(j)
        if entries:
            for row_idx, coeff in entries:
                raw_met = str(mets[row_idx])
                met_obj = met_objs.get(raw_met)
                if met_obj is None:
                    fallback_id, fallback_comp = _transform_matlab_met_id(raw_met)
                    if fallback_id not in existing_met_ids:
                        fallback_met = cobra.Metabolite(fallback_id, compartment=fallback_comp)
                        model.add_metabolites([fallback_met])
                        existing_met_ids.add(fallback_id)
                    met_obj = model.metabolites.get_by_id(fallback_id)
                rxn.add_metabolites({met_obj: float(coeff)})

        if j < len(cvec):
            rxn.objective_coefficient = float(cvec[j])

    if genes is not None:
        for gene_name in genes:
            gene_id = _normalize_gene_id(str(gene_name))
            if gene_id:
                model.genes.get_by_id(gene_id)

    if gr_rules is not None:
        for rxn, rule in zip(model.reactions, gr_rules):
            if rule is None:
                continue
            rule = _normalize_gpr(rule)
            if rule:
                rxn.gene_reaction_rule = rule

    return model


def prepare_mcf7() -> None:
    """Prepare the MCF7 model from the checked-in SBML file."""
    print("Preparing MCF7 model")

    output_model_path = MODELS_DIR / "MCF7_Zielinski_2017_core.sbml"
    if output_model_path.exists():
        model = read_sbml_model(str(output_model_path))
        print(f"Using existing MCF7 model: {output_model_path}")
        print(f"- {model.id}: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
        return

    if loadmat is None or csc_matrix is None:
        raise RuntimeError(
            "SciPy is required to regenerate the MCF7 model from the MATLAB archive. "
            "The repository already ships the prepared SBML file at models/MCF7_Zielinski_2017_core.sbml."
        )

    matfile = MODELS_DIR / "Zielinski_et_al_2017" / "code" / "modelsAll.mat"
    label_file = MODELS_DIR / "Zielinski_et_al_2017" / "cell_line_labels.txt"
    if not matfile.exists():
        raise FileNotFoundError(f"Missing MATLAB model archive: {matfile}")

    data = loadmat(str(matfile), struct_as_record=False, squeeze_me=True)
    models = data["modelsAll"]

    with label_file.open() as handle:
        cell_line_labels = [line.rstrip() for line in handle]

    assert len(models) == len(cell_line_labels)

    cobra_model_list = []
    for index, model_struct in enumerate(models):
        cobra_model = _cobramat2cobrapy(model_struct)
        cobra_model.id = cell_line_labels[index]
        cobra_model_list.append(cobra_model)

    mcf7_model = next(model for model in cobra_model_list if model.id == "MCF7")

    for reaction in mcf7_model.reactions:
        if reaction.objective_coefficient != 0:
            reaction.objective_coefficient = 0
    mcf7_model.reactions.biomass_NCI60.objective_coefficient = 1

    solution = mcf7_model.optimize()
    print(f"MCF7 objective: {solution.objective_value}")

    write_sbml_model(mcf7_model, str(output_model_path))
    print(f"Wrote {output_model_path}")


def main() -> None:
    prepare_clostridium()
    prepare_methanosarcina()
    prepare_iml1515()
    trim_models()
    prepare_mcf7()


if __name__ == "__main__":
    main()
