# Rat_TTR_PFAS_QSAR_Docking
Repository associated with the manuscript: "New Docking, Molecular Dynamics and QSAR Models to Predict Disruption of Human and Rat Transthyretin Function by Per- and Polyfluoroalkyl Substances (PFAS)". This repo contains the analysis code, input data, and files used to support the study.

## Repository map
- `R/Bayesian_regression_TTR_inVivoPODs.qmd` - Quarto analysis for linear and Bayesian regression of rat TTR binding energy vs in vivo PODs; reads input data, transforms units, and generates plots.
- `R/dose_to_serum_pfas.Rmd` - R Markdown workflow to convert administered dose to serum metrics using TK parameters and study-specific exposure data.
- `R/dose_serum_functions.R` - Helper functions used by `R/dose_to_serum_pfas.Rmd` for PK/TK calculations and simulations.
- `data_input/Thyroid_dataset_1.22.2025.xlsx` - Source dataset with BMDs and binding energies (see "Condensed" and "binding" sheets).
- `data_input/T4_model_averaging_results_molar_conversions correct file.csv` - Updated T4 model averaging POD conversions used to override uM BMDs.
- `figures/ - Figure output produced by R and Python used in the manuscript.

## Reproducibility
- R version 4.5.1 was used for these analysis, however a recent R 4.x release is expected to work.
- R packages used in `R/Bayesian_regression_TTR_inVivoPODs.qmd`: `tidyverse`, `cols4all`, `ggpmisc`, `ggpubr`, `rstanarm`, `bayesplot`, `brms`, `dplyr`, `broom`, `purrr`, `scales`, `readxl`, `readr`.
- Stan toolchain: `rstanarm`/`brms` require a C++ toolchain (e.g., Rtools on Windows) and Stan configured for your system.
 - R package versions: writexl_1.5.4, scales_1.4.0, broom_1.0.11, brms_2.23.0, bayesplot_1.15.0, rstanarm_2.32.2, Rcpp_1.1.0, ggpubr_0.6.2, ggpmisc_0.6.3, ggpp_0.5.9, cols4all_0.10, lubridate_1.9.4, forcats_1.0.1, stringr_1.6.0, dplyr_1.1.4, purrr_1.2.1, readr_2.1.6, tidyr_1.3.2, tibble_3.3.0, ggplot2_4.0.1, tidyverse_2.0.0 

