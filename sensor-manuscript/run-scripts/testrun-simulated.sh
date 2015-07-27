#!/bin/bash
#module load R/3.2.0
sensor=${1:-'aviris.ng'}
ngibbs=${2:-'500'}
qsub -V -v N=1.068680,Cab=28.006041,Car=6.323790,Cw=0.003950,Cm=0.000455,sensor=$sensor,ngibbs=$ngibbs,runname=N1.q1.$sensor -N "testrun.$sensor" submit.simulated.qsub
