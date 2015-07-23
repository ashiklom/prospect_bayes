#!/bin/bash
# Arguments: sensor ngibbs
sensor=${1:-identity}
ngibbs=${2:-50000}

while read l; do
    qsub -V -v spectra=$l,ngibbs=$ngibbs,sensor=$sensor -N "$l $sensor inversion" submit.fft.qsub
done < fft-speclist.dat
