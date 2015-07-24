##' Load traits into table
library(reshape2)
library(data.table)

# Read command line arguments
arg <- commandArgs(trailingOnly=TRUE)
PATH.data <- arg[1]
if(is.na(PATH.data)) PATH.data <- "/mnt/dropbox/NASA_TE_PEcAN-RTM_Project/Data"
PATH.spectra <- arg[2] 
if(is.na(PATH.spectra)) PATH.spectra <-"/mnt/dropbox/FFT_spectra/FFT_spectra_unique.csv"
PATH.speciesinfo <- arg[3]
if(is.na(PATH.speciesinfo)) PATH.speciesinfo <- "/mnt/dropbox/NASA_TE_PEcAn-RTM_Project/Data/FFT_species_info_csv.csv"

### Setup paths
PATH.data.FFT = file.path(PATH.data, "NASA_FFT_Project")
PATH.data.FFT.d15N = file.path(PATH.data.FFT, "NASA_FFT_d15N_ANALYZED_DATA_UPDATED_4R.csv")
PATH.data.FFT.lignin = file.path(PATH.data.FFT, "NASA_FFT_FIB_LIG_CELL_RESULTS_FINAL_4R.csv")
PATH.data.FFT.CN = file.path(PATH.data.FFT, "NASA_FFT_Project_CN_Data_4R.csv")
PATH.data.FFT.SLA_LMA = file.path(PATH.data.FFT, "NASA_FFT_SLA_LMA_Data_v2_4R_updated_new.csv")
PATH.results <- "data/fft.dat.RData"

### Read data
FFT.d15N <- read.csv(PATH.data.FFT.d15N, header=TRUE, stringsAsFactors = FALSE)
FFT.lignin <- read.csv(PATH.data.FFT.lignin, header=TRUE, stringsAsFactors = FALSE)
FFT.CN <- read.csv(PATH.data.FFT.CN, header=TRUE, stringsAsFactors = FALSE)
FFT.SLA <- read.csv(PATH.data.FFT.SLA_LMA, header=TRUE, stringsAsFactors = FALSE)

### Remove bad data
FFT.d15N <- subset(FFT.d15N, COMMENTS == "")
FFT.lignin <- subset(FFT.lignin, COMMENTS == "")
FFT.SLA$EWT_g_cm2[FFT.SLA$EWT_g_cm2 < 0] <- NA
FFT.SLA$LMA_g_DW_cm2[FFT.SLA$LMA_g_DW_cm2 < 0] <- NA

### Subset necessary columns
mergeby.caps <- c("SAMPLE_NAME",
                  "SAMPLE_YEAR")
mergeby.lower <- c("Sample_Name",
                   "Sample_Year")
FFT.d15N.sub <- FFT.d15N[,c(mergeby.caps,
                            "SAMPLE_dN15")]
FFT.lignin.sub <- FFT.lignin[,c(mergeby.caps,
                                "ADF_PERC_DW",
                                "ADL_PERC_DW",
                                "CELL_PERC_DW")]
FFT.CN.sub <- FFT.CN[, c(mergeby.lower,
                         "Perc_N",
                         "Perc_C",
                         "CNRatio")]
FFT.SLA.sub <- FFT.SLA[, c(mergeby.lower,
                           "EWT_g_cm2",
                           "LMA_g_DW_cm2")]

### Merge into large data file
FFT.p1 <- merge(FFT.d15N.sub, FFT.lignin.sub, by=mergeby.caps, all=TRUE)
FFT.p2 <- merge(FFT.CN.sub, FFT.SLA.sub, by=mergeby.lower, all=TRUE)
FFT.all <- merge(FFT.p2, FFT.p1, by.x=mergeby.lower, by.y=mergeby.caps, all=TRUE)
FFT.all <- data.table(FFT.all)

# Load inversion results ("fft.dat")
load(PATH.results)

## Merge data tables
info.orig <- fread(PATH.spectra, header=TRUE, select = c(1:9,12,13))
setkey(info.orig, Spectra)
setkey(fft.dat, spectra)
results <- fft.dat[info.orig]

## Average duplicates in chemistry data
FFT.nodup <- FFT.all[, lapply(.SD, mean, na.rm=TRUE), by=c("Sample_Name", "Sample_Year")]

## Merge results with traits
fftdat <- merge(x=results, y=FFT.nodup, by=mergeby.lower, all=TRUE)

## Sort out missing values
sample.name.regex <- "^([A-Za-z]+)([0-9]+[B]*)_([A-Za-z]+)_([BMTS]+)([ANO23]{0,1})[SAMP2]*"
fftdat[is.na(Species), Site := gsub(sample.name.regex, "\\1", Sample_Name)]
fftdat[is.na(Species), Plot := gsub(sample.name.regex, "\\1\\2", Sample_Name)]
fftdat[is.na(Species), Age := gsub(sample.name.regex, "\\5", Sample_Name)]
fftdat[is.na(Species), Height := gsub(sample.name.regex, "\\4", Sample_Name)]
fftdat[is.na(Species), Species := gsub(sample.name.regex, "\\3", Sample_Name)]

## Load species information
species.info <- fread(PATH.speciesinfo, header=TRUE)
species.info[,Scientific := sprintf("%s %s", Genus, Species)]
species.info[, Species := NULL]
setnames(fftdat, "Species", "Label")
setkey(fftdat, "Label")
setkey(species.info, "Label")
uk <- unique(c(fftdat[,Label], species.info[,Label]))
fft.full <- fftdat[species.info[J(uk)]]
save(fft.full, file="data/FFT_full.RData")
