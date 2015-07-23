#!/bin/bash
#module load R/3.2.0
testspec=${1:-'AK01_ACRU_B_LC_REFL_2009'}
qsub -V -v spectra=$testspec,ngibbs=800,sensor=identity -N "testrun" submit.fft.qsub
