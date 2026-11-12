# dfba_sbml_models

This repository contains the genome-scale metabolic models and preprocessing scripts used in the **PhysiCelldFBA** manuscript. It archives the exact model versions employed in the simulation case studies together with the scripts required to retrieve, standardize, curate, and validate them from their original published reconstructions.

## Repository structure

```
.
├── models/                               Original and curated metabolic models
├── modules/                              Python helper functions
├── notebooks/                            Preparation notebooks for each model
│   ├── C_beijerinckii_model_iCB925.ipynb
│   ├── M_barkeri_iMG746_preparation.ipynb
│   └── load_Zielinski_2017_core_human_models.ipynb
├── prepare_models.py                     Orchestrates the full curation workflow
├── retrieve_models.sh                     Download publicly available models
├── RPMI_1640_media.json                   RPMI-1640 medium composition
├── requirements.txt
├── LICENSE
└── README.md
```

## Metabolic models

The repository contains the genome-scale metabolic models used in the different PhysiCelldFBA simulation examples.

### *Escherichia coli* core metabolic model

**Reference**

Orth JD, Fleming RMT, Palsson BØ.
*Reconstruction and use of microbial metabolic networks: the core Escherichia coli metabolic model as an educational guide.*
EcoSal Plus. 2010.

**Source**

BiGG Models

https://bigg.ucsd.edu/models/e_coli_core

---

### iML1515 (*Escherichia coli*)

**Reference**

Monk JM, Lloyd CJ, Brunk E, et al.
*iML1515, a knowledgebase that computes Escherichia coli traits.*
Nature Biotechnology. 2017;35(10):904–908.

**Source**

BiGG Models

https://bigg.ucsd.edu/models/iML1515

---

### iCB925 (*Clostridium beijerinckii* NCIMB 8052)

**Reference**

Milne CB, Eddy JA, Raju R, Ardekani S, Kim PJ, Senger RS, Jin YS, Blaschek HP, Price ND.
*Metabolic network reconstruction and genome-scale model of butanol-producing strain Clostridium beijerinckii NCIMB 8052.*
BMC Systems Biology. 2011;5:130.

**Source**

BioModels

https://www.ebi.ac.uk/biomodels/MODEL1507180034

Repository files

- `models/C_iCB925_beijerinckii.xml` (original reconstruction)
- `models/C_beijerinckii_model_iCB925_trimmed.xml` (curated model used in the manuscript)

---

### iMG746 (*Methanosarcina barkeri* Fusaro)

**Reference**

Gonnerman MC, Benedict MN, Feist AM, Metcalf WW, Price ND.
*Genomically and biochemically accurate metabolic reconstruction of Methanosarcina barkeri Fusaro, iMG746.*
Biotechnology Journal. 2013;8(9):1070–1079.

**Source**

Supporting Information S2 of the original publication.

Repository files

- `models/M_barkeri_model_iMG746.xml`
- `models/M_barkeri_model_iMG746_trimmed.xml`

Note on the trimmed iCB925 and iMG746 files (changed in this version; previous files: commit 3d95bdf)

The trimmed models of the two species used in the cross-feeding example had every remaining
reaction bounded to its flux range in the medium used to prepare them (flux caps). For iMG746 these
caps restrict the model in the cross-feeding medium (518 internal bounds changed, 314 of them upper
bounds cut from 1000 to below 1, optimum of the preparation medium 15.7% lower than in the full
model). In the files of this version the reactions removed by the trimming (blocked in the
preparation medium) and the stoichiometry, gene associations and objective are unchanged, and every
internal reaction has the bounds of the full model (`models/M_barkeri_model_iMG746.xml` and the
full iCB925 reconstruction) instead of the flux caps (485 bounds changed in iCB925, 518 in iMG746).
Three further differences: the upper bound of the biomass reaction (0.35 per h in iCB925 and 0.05 per h
in iMG746, set in the preparation) is removed, so growth is limited by the uptake only; the ATPM bounds
are those of the full model with the lower bound kept at 0 ((0, 300) to (0, 1000) in iCB925 and (0, 20)
to (0, 2) in iMG746); and in iMG746 the uptake of HCO3-, SO3, N2 and H+ is closed (exchange lower
bound 0) because no cell used it in the simulations of the cross-feeding example. `prepare_models.py`
generates these files (`model_trimming(..., cap_internal_bounds=False)` is the default). These are the models used in the revised manuscript; the files and the
simulations that use them are archived with it (see the manuscript for the DOI).

sha256 (first 12 characters): `C_beijerinckii_model_iCB925_trimmed.xml` 76c8741fbf4c,
`M_barkeri_model_iMG746_trimmed.xml` 93465969ef09.

---

### MCF7 genome-scale metabolic model

**Reference**

Zielinski DC, Jamshidi N, Corbett AJ, Bordbar A, Thomas A, Palsson BØ.
*Systems biology analysis of drivers underlying hallmarks of cancer cell metabolism.*
Scientific Reports. 2017;7:41241.

**Source**

Supplementary Information 3 of the original publication.

Repository file

- `models/MCF7_Zielinski_2017_core.sbml`

The repository also includes the supplementary material from the original publication containing the complete collection of reconstructed cancer cell-line metabolic models.

---

## Downloading the models

The publicly available models can be downloaded automatically by running

```bash
bash retrieve_models.sh
```

This script downloads:

- **iCB925** from BioModels.
- **iML1515** from the BiGG database.
- **E. coli core** from the BiGG database.

The **iMG746** reconstruction is distributed as supplementary material of Gonnerman *et al.* (2013) and must be downloaded manually following the instructions printed by the script.

---

## Model curation

The models distributed in this repository correspond to the curated versions used in the PhysiCelldFBA manuscript.

The preprocessing workflow performs several model-specific operations, including:

- Standardization of reaction and metabolite identifiers.
- Standardization of exchange reaction naming conventions.
- Correction of SBML compatibility issues.
- Configuration of simulation-specific environmental conditions.
- Harmonization of exchange reaction bounds across interacting models.
- Removal of blocked reactions.
- Removal of orphan (gap) metabolites.
- Validation that the curated models preserve the original optimal biomass objective value.

The resulting curated models are written to the `models/` directory.

---

## Reproducing the preprocessing

Create a Python environment and install the required dependencies:

```bash
python -m venv .venv
source .venv/bin/activate

pip install -r requirements.txt
```

Then execute

```bash
python prepare_models.py
```

Model-specific preparation steps are documented in the accompanying notebooks under the notebooks directory:

- `notebooks/C_beijerinckii_model_iCB925.ipynb`
- `notebooks/M_barkeri_iMG746_preparation.ipynb`
- `notebooks/load_Zielinski_2017_core_human_models.ipynb`

These notebooks describe the preprocessing performed to generate the curated models used in the PhysiCelldFBA simulation examples.

---

## Attribution

The original genome-scale metabolic models were developed by their respective authors and are redistributed or referenced according to their original licenses and repositories.

This repository provides the curated versions used in the PhysiCelldFBA publication together with the preprocessing scripts required to reproduce the curation workflow. The original authors retain all rights to their respective metabolic reconstructions.

---

## Citation

If you use the curated models or preprocessing workflow provided in this repository, please cite:

- The original publications describing each metabolic reconstruction.
- The **PhysiCelldFBA** publication.
- The Zenodo archive associated with this repository.

---

## License

This repository is distributed under the terms specified in the `LICENSE` file.
