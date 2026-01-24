# Rat_TTR_PFAS_QSAR_Docking
Repository associated with the manuscript: "New Docking, Molecular Dynamics and QSAR Models to Predict Disruption of Human and Rat Transthyretin Function by Per- and Polyfluoroalkyl Substances (PFAS)". This repo contains the analysis code, input data, figures, and manuscript files used to support the study.

## Repository map
- `Main_document_TTR_PFAS.docx` - Manuscript draft.
- `R/Bayesian_regression_TTR_inVivoPODs.qmd` - Quarto analysis for linear and Bayesian regression of rat TTR binding energy vs in vivo PODs; reads input data, transforms units, and generates plots.
- `data_input/Thyroid_dataset_1.22.2025.xlsx` - Source dataset with BMDs and binding energies (see "Condensed" and "binding" sheets).
- `data_input/T4_model_averaging_results_molar_conversions correct file.csv` - Updated T4 model averaging POD conversions used to override uM BMDs.
- `figures/figure9.jpg` - Figure output used in the manuscript.
- `Rat_TTR_PFAS_QSAR_Docking.Rproj` - RStudio project file.
- `LICENSE` - Repository license.

## Reproducibility
- R version: not recorded in the repository; capture your exact version with `sessionInfo()` for reporting. A recent R 4.x release is expected to work.
- R packages used in `R/Bayesian_regression_TTR_inVivoPODs.qmd`: `tidyverse`, `cols4all`, `ggpmisc`, `ggpubr`, `rstanarm`, `bayesplot`, `brms`, `dplyr`, `broom`, `purrr`, `scales`, `readxl`, `readr`.
- Stan toolchain: `rstanarm`/`brms` require a C++ toolchain (e.g., Rtools on Windows) and Stan configured for your system.
- Working directory: run analyses from the repository root so relative paths in `R/Bayesian_regression_TTR_inVivoPODs.qmd` resolve to `data_input/`.
