## Run of Bayesian inversion of PROSPECT

args <- commandArgs(trailingOnly=TRUE)
jrsd <- as.numeric(args[1])
ngibbs <- as.numeric(args[2])
runid <- args[3]
folder <- args[4]
filename <- sprintf("run_results/%s/ALL_%g_hre_%s.dat", folder, jrsd, runid)
dir.create(sprintf("run_results/%s", folder), showWarnings = FALSE)

source("inv_bayes.R")
source("specdataproc.R")

smat <- specmatrix("All")

pinvbayes(smat, ngibbs=ngibbs, JumpRSD=jrsd, fname=filename,
          inits='random',
          ar.step=100,
          notifications=FALSE)
