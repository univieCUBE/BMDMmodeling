THis is the folder containing custom scripts used for the purposes of the manuscript, listed below in the chronological order of the analysis:

1. `TPM_prepINIT.R` - code for transforming gene count values into TPM values, to be used for contextualizing metabolic models with MOUSE GEM as template
2. `tINIT2extraction_Mouse_GEM.m` - code for contextualizing metabolic models using tINIT algorithm
3. `nontarg.R` - code for calculating uptake and secretion rates for metabolites measured in the supernatant in an untargeted fashion, for constraining purposes
4. `ExchangeMetabolomics.m` - code for incorporating calculated uptake and secretion rates into the models. Afterwards, fluxes are simulated with FBA, FVA and random sampling.
5. `Main.R` - code for all downstream purposes post-flux prediction: vizualization of the figures as shown in the manuscript, as well as statistics for sampling results. Lastly, code for analyzing validation of intracellular metabolomics is shown.
