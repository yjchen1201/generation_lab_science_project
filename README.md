# generation_lab_science_project

This repository contains the code and data for the DNA methylation project. The goal of this project is to ensure consistency across DNA methylation data from various public datasets and an internal dataset.

## Project Overview

The project involves the following steps:
1. Load and preprocess data.
2. Perform quality control.
3. Normalize the data.
4. Correct batch effects.
5. Further quality checks.
6. Save corrected data.

## Files

- `Generation_Lab_proj.R`: The main R script containing all the steps.
- `public_beta_values_corrected`: Corrected beta value csv file for each public sample.
- `internal_beta_values_corrected`: Corrected beta value csv file for each internal sample.
- `Generation_lab_proj_sessionInfo.txt`: Information about the R session, including package versions.

## Usage

To reproduce the results, run the `Generation_Lab_proj.R` script in R.

```R
source("Generation_Lab_proj.R")
