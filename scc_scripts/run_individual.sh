#!/bin/bash -l

## Perform an inversion. Parameters are as follows:
## [1] Spectra
## [2] Number of iterations
## [3] Folder name

module load R/R-3.1.1

SPECARG=$1
NGIBBSARG=${2:-1000}
FOLDARG=${3:-testfolder}

qsub -v SPECTRA=$SPECARG,NGIBBS=$NGIBBSARG,FOLDER=$FOLDARG,RUN=01 -N "pbi $SPECARG $FOLDARG $RUN" submit_prosp_nimble.qsub
qsub -v SPECTRA=$SPECARG,NGIBBS=$NGIBBSARG,FOLDER=$FOLDARG,RUN=02 -N "pbi $SPECARG $FOLDARG $RUN" submit_prosp_nimble.qsub
qsub -v SPECTRA=$SPECARG,NGIBBS=$NGIBBSARG,FOLDER=$FOLDARG,RUN=03 -N "pbi $SPECARG $FOLDARG $RUN" submit_prosp_nimble.qsub
qsub -v SPECTRA=$SPECARG,NGIBBS=$NGIBBSARG,FOLDER=$FOLDARG,RUN=04 -N "pbi $SPECARG $FOLDARG $RUN" submit_prosp_nimble.qsub
qsub -v SPECTRA=$SPECARG,NGIBBS=$NGIBBSARG,FOLDER=$FOLDARG,RUN=05 -N "pbi $SPECARG $FOLDARG $RUN" submit_prosp_nimble.qsub

