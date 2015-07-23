#!/bin/bash
#module load R/3.2.0
testspec=${1:-'AK01_ACRU_B_LC_REFL_2009'}
sensor=${2:-'aviris.ng'}
qsub -V -v spectra=$testspec,ngibbs=800,sensor=$sensor -N "testrun" submit.fft.qsub
