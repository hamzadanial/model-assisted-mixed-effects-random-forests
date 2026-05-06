# Model-Assisted Small Area Estimation Using Mixed-Effects Random Forests

This repository contains the R code and supporting datasets used for the analyses in the manuscript:

**“Model-Assisted Small Area Estimation Using Mixed-Effects Random Forests”**

The repository includes code for:

1. simulation data generation,
2. simulation study analysis,
3. PDHS data analysis,
4. bootstrap MSE analysis for MERF and MA-MERF under the simulation study,
5. bootstrap MSE analysis for MERF and MA-MERF under the PDHS application.

---

## Repository contents

### R scripts
- `01_simulation_data_generation.R`  
  Generates the simulated finite population used in the simulation study.

- `02_simulation_study_analysis.R`  
  Performs the main simulation study analysis and produces the corresponding summary results.

- `03_pdhs_data_analysis.R`  
  Conducts the empirical analysis using the refined PDHS dataset.

- `04_bootstrap_mse_simulation_merf_mamerf.R`  
  Computes bootstrap-based MSE measures for MERF and MA-MERF under the simulation study.

- `05_bootstrap_mse_pdhs_merf_mamerf.R`  
  Computes bootstrap-based MSE measures for MERF and MA-MERF under the PDHS-based bootstrap study.

### Data files
- `simulated_data.csv`  
  Simulated dataset used in the simulation study.


Due to data-use restrictions, the PDHS/DHS dataset is not included in this repository. Researchers who wish to reproduce the analysis may obtain the data directly from The DHS Program and place the required files locally before running the scripts. 

For questions about the code structure or data preparation steps, they may contact the corresponding author.
---

## Software used

The analyses were carried out using:

- **R** version **4.4.3**
- **RStudio Server** version **2024.12.1 Build 563**

---

## Main R packages

The main R packages used in this repository include:

- `data.table`
- `emdi`
- `SAEforest`
- `lme4`
- `ranger`
- `parallel`
- `pbapply`
- `reshape2`

Additional package requirements, if any, are specified within the corresponding scripts.

---

## Recommended order of execution

To reproduce the analyses, the scripts may be run in the following order:

1. `01_simulation_data_generation.R`
2. `02_simulation_study_analysis.R`
3. `04_bootstrap_mse_simulation_merf_mamerf.R`
4. `03_pdhs_data_analysis.R`
5. `05_bootstrap_mse_pdhs_merf_mamerf.R`

---

## Notes

- The simulation study and the PDHS-based analysis are provided separately for clarity and reproducibility.
- The bootstrap scripts are written specifically for comparing **MERF** and **MA-MERF** under both simulation and PDHS settings.
- File names and script order were organized to make the workflow easy to follow.

---

## Correspondence

For questions regarding the code or analysis workflow, please contact the authors of the manuscript.
