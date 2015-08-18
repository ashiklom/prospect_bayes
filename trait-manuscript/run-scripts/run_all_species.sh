#! /bin/bash
## Read in SE species and pass to 'species_run.sh'
##
## Arguments to species_run.sh:
## [1] Spectra
## [2] Number of iterations
## [3] Folder

NG=$1
FD=$2

while read l; do
        echo $l
        run_species.sh $l $NG $FD
done < ../data/species_list_FFT.txt
