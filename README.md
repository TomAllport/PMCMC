# Particle Markov chain Monte Carlo methods for joint state and parameter estimation (Year-Long Project)

This repository contains two R files to accompany my MA40249 year-long project.

## Required Packages

The following packages are required to run the files:

```{r}
library('abind')
library('MASS')
library('tictoc')
library('grid')
```

## Usage

The code provided reproduces the analysis and plots presented in my project report. For reproducibility, the seed has been set to the same value in both files. Some of the sections of ```PMH.R``` take considerable time to run so is recommended to only run these when necessary. These are lines 111-122 (to calculate the time taken for the PMMH to run for different numbers of particles) and lines 124-138 (to calculate the acceptance for the PMMH for different numbers of particles and proposal standard deviations).
