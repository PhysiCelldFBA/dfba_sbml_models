# metabolic_models
A repository for fixing and tweaking the metabolic models used in the use cases of the paper


Model iCB925
- Metabolic network reconstruction and genome-scale model of butanol-producing strain Clostridium beijerinckii NCIMB 8052.
- Milne CB, Eddy JA, Raju R, Ardekani S, Kim PJ, Senger RS, Jin YS, Blaschek HP, Price ND.
- BMC Syst Biol 2011; 5: 130

Model retrived from BioModels: https://www.ebi.ac.uk/biomodels/MODEL1507180034


Model iAF692
- Modeling methanogenesis with a genome-scale metabolic reconstruction of Methanosarcina barkeri
- Adam M Feist 1, Johannes C M Scholten, Bernhard Ø Palsson, Fred J Brockman, Trey Ideker
- iMol Syst Biol. 2006:2:2006.0004.

Model retrived from BiGG database: http://bigg.ucsd.edu/models/iAF692


## Retriving models
To retrive models run `bash retrive_models.sh` 

## Trimming models
To prepare the models to be used in the the methanogenesis from CO2 & H2 crossfeeding use case we applied the following steps:
1- Fixing reactions and models ID (when needed)
2- Set boundory conditions in the two models to perform the methanogenesis from CO2 & H2 crossfeeding use case. 
3- Models are trimmed by removing blocked reactions and gap metabolties.
4- Trimmed models are tested to to have an objective biomass value equalt to the one in the original model
5- Models trimeed and fixed are stored in the models folder  


  
