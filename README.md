# Arbo_riskmaps_public

Repository contains data, scripts, and key outputs from "The overlapping global distribution of dengue, chikungunya, Zika, and yellow fever".

## System Requirements

This code was developed and tested using:

- **R Version:** 4.3.0
- **Operating System:** Windows 10 x64 (build 19045)

### Software Dependencies

Ensure the following R packages are installed: `pacman (v0.5.1)`, `doParallel` (v1.0.17), `foreach` (v1.5.2), `pdp` (v0.8.1), `matrixStats` (v1.3.0), `Metrics` (v0.1.4), `cutpointr` (v1.1.2), `boot` (v1.3-28.1), `blockCV` (v3.1-3), `rsample` (v1.2.0), `randomForest` (v4.7-1.1), `terra` (v1.7-29), `rnaturalearth` (v0.3.4), `raster` (v3.6-23), `sp` (v1.6-0), `exactextractr` (v0.9.1), `sf` (v1.0-14), `tidyterra` (v0.5.2), `countrycode` (v1.5.0), `dplyr` (v1.1.2), `data.table` (v1.14.8). 

### Non-Standard Hardware

No special hardware requirements beyond a standard computer with sufficient RAM to handle the scale of analysis.

### Tested Versions

- **R Version:** 4.3.0
- **Operating System:** Windows 10 x64 (build 19045)

## Installation Guide

### Instructions

1. **Install R:** Download and install R from [CRAN](https://cran.r-project.org/).
2. **Install Required Packages:** Run the following R script to install the required packages:
   ```r
   if(!require("pacman")) install.packages("pacman")
   # list of the required packages: 
   pkgs = 
      c(# data manipulation
        "tidyr", "data.table","dplyr", "countrycode", "conflicted", "stringi", "lubridate",
        # visualizations
        "ggplot2", "patchwork", "ggpubr", "viridis", "ggtext", "glue", "extrafont", "tidyterra",
        # spatial data
        "sf", "mapview", "exactextractr", "raster", "rnaturalearth", "terra", 
        # modelling 
        "randomForest", "rsample", "blockCV", 
        "boot", "cutpointr", 
        "Metrics", "pROC", "ROCR", "matrixStats", "pdp",
        # programming 
        "here", "tictoc", "doParallel"
      ) 
   pacman::p_load(pkgs, character.only = T)

### Typical Install Time
The installation of R and the required packages typically takes less than 30 minutes on a standard desktop computer.

## Instructions for Use

### Repository Structure

The repository is organized into the following main folders:

- **`data/`**: Contains datasets used in the analysis, subdivided into:
  - **`admin_rasters/`**: Rasters and shapefiles for global administrative boundaries.
  - **`covariate_rasters/`**: Rasters for global dengue temperature suitability, population density, and a mask layer for yellow fever transmission. See Methods and Supplementary Information for data sources.
  - **`intermediate_datasets/`**: Cleaned and standardized occurrence datasets for viral and arboviral diseases. See cited sources for original datasets. This folder also contains several intermediate datasets that are used at various stages of the analysis. 

- **`functions/`**: Contains custom R code used in various stages of the analysis. Load these functions into your R script before use.

- **`script/`**: Contains all analysis code and scripts for generating figures. Key scripts include:
  - `00_setup.R`: Sets up R packages and functions.
  - `01_get_data_arbo_model.R`, `01_get_data_surv_model.R`: Prepares covariate and disease occurrence data.
  - `02_surv_model_fitting.R`, `02_arbo_model_fitting.R`: Fits the surveillance capability and arbovirus models.
  - `04_masking_and_binary_maps.R`: Post-processes predicted maps of disease suitability.
  - `05_pop_at_risk.R`: Calculates global and regional population at risk of arboviral diseases.
  - `06_visualisations.R`, `07_supplementary_figs.R`: Generates main and supplementary figures.

- **`outputs/`**: Contains outputs obtained from the analyses, subdivided into:
  - **`Cross_validation/`**: Contains out-of-bag predictions, model diagnostics, and sensitivity analysis results.
  - **`Tables/`**: Contains tables obtained from the analysis.
  - **`Figures/`**: Contains all main and supplementary figures in the manuscript.
  - **`Rasters/`**: Contains all rasters obtained from analysis and used to create figures in the manuscript.

For additional details, refer to the comments within individual scripts and any supplementary documentation provided.

