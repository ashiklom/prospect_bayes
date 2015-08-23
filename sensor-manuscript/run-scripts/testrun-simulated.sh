#!/bin/bash
#module load R/3.2.0
sensor=${1:-'chris.proba'}
ngibbs=${2:-'5000'}
qsub -V -v index=1,sensor=$sensor,ngibbs=$ngibbs,runname=testrun.1.$sensor -N "testrun.$sensor" submit.simulated.qsub
