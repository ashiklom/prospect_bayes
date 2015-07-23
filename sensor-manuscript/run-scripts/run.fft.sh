#!/bin/bash -l
# Arguments: ngibbs
module load R/3.2.0
sensor=$(1:-identity)
ngibbs=$(2:-50000)

while read l; do
    qsub -v spectra=$l,ngibbs=$ngibbs,sensor=$sensor -N "$l $sensor inversion" submit.fft.qsub
done < fft-speclist.dat
