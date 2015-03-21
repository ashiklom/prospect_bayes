#!/bin/bash

## Script for test inversion of PROSPECT model.
## Arguments to "invert_individual.R" are as follows:
##    [1] Spectra
##    [2] Number of MCMC steps 
##    [3] Folder for result storage
##    [4] Run ID

touch ../run_results/testfolder/empty
rm ../run_results/testfolder/*
Rscript invert_individual.R AK01_ACRU_B_LC_REFL_2009 100 test_folder test1
