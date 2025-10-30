## 16S Sequencing Analysis Pipeline for Low Biomass Samples

This repository contains an analysis pipeline for 16S sequencing data developed for low biomass microbiome samples.  
All code and methods were developed and are maintained by Jessica Whelan (jessica-whelan).
The pipeline is modular and under active development. 
---

## Repository Contents

This repository contains the R code and analysis pipeline for 16S microbiome data analysis used in the **bmd-project**. Below is a description of the main scripts and their intended purpose:

- **`16S-sequencing.R`** – The main analysis script for the project.  
  - Performs preprocessing of raw 16S sequencing data, including quality filtering, decontamination, and normalization.  
  - Computes alpha and beta diversity metrics to assess microbial diversity within and between samples.  
  - Implements ML-based differential abundance analyses, allowing identification of taxa associated with experimental conditions.  
  - Serves as the central workflow for the project; additional functions and plots may be integrated as development progresses.

- **`krona-plot-function.R`** – Interactive visualization module (in development).  
  - Generates Krona HTML plots for hierarchical exploration of taxonomic abundance.  
  - Designed to provide an intuitive, interactive view of microbial composition for presentations or exploratory analysis.  
  - Once finalized, outputs will be  incorporated into the main analysis script.

- **`testing-code.R`** – Development and testing script.  
  - Contains code for various visualization options and exploratory analyses.  
  - Includes **ANCOM-BC2** analysis for sensitivity testing of differential abundance results.
  - Includes intra-sample variability analyses 
  - Functions here are experimental and may be merged into the main workflow once validated.

Notes:
- All scripts are written in **R** and rely on packages such as `phyloseq`, `mixOmics`, and other common microbiome analysis tools. Please see CITATION.md for references to R packages used in this pipeline that request citation

---

## Authorship and Usage

This pipeline was developed and is maintained by Jessica Whelan (jessica-whelan). It serves as the upstream source for the analysis code used in the bmd-project repo.  
A corresponding fork is maintained under the `bmd-project` organization to support project-specific development and collaboration.

Please acknowledge or cite this repo when using or adapting the code.
Any use, adaptation, redistribution, or incorporation of the code requires permission from the author and must include appropriate citation or acknowledgement.

© 2025 Jessica Whelan. All rights reserved.
