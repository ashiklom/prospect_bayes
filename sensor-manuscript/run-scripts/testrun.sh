#!/bin/bash
#module load R/3.2.0
sensor=${1:-'aviris.ng'}
testspec=${2:-'AK01_ACRU_B_LC_REFL_2009'}
qsub -V -v spectra=$testspec,ngibbs=800,sensor=$sensor -N "testrun" submit.fft.qsub
