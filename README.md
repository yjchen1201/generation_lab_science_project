# generation_lab_science_project

This repository contains the code and data for the DNA methylation project. The goal of this project is to ensure consistency across DNA methylation data from various public datasets and an internal dataset. The key outputs are corrected beta values for each sample.

## Project Overview

The project involves the following steps:
1. Load and preprocess data.
2. Perform quality control.
3. Normalize the data.
4. Correct batch effects.
5. Further quality checks.
6. Save corrected data.

## Files

- `Generation_lab_proj_sessionInfo.txt`: Information about the R session, including package versions.
- `Generation_Lab_proj.R`: The main R script containing all the steps.
- `public_beta_values_corrected`: Corrected beta value csv files for each public sample. https://drive.google.com/drive/folders/1dkrOGaLcCbnz5yW-ObJlm-kest9WnW_K?usp=sharing
- `internal_beta_values_corrected`: Corrected beta value csv files for each internal sample. https://drive.google.com/drive/folders/1PBCAKRcHvGHgxHGr8Yd5cC3f3gOOt36U?usp=sharing

