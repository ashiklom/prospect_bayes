##' Load traits into table
library(reshape2)
library(data.table)

### Setup paths
PATH.data = "~/Documents/Dropbox/NASA_TE_PEcAn-RTM_Project/Data"
PATH.data.FFT = file.path(PATH.data, "NASA_FFT_Project")
PATH.data.FFT.d15N = file.path(PATH.data.FFT, "NASA_FFT_d15N_ANALYZED_DATA_UPDATED_4R.csv")
PATH.data.FFT.lignin = file.path(PATH.data.FFT, "NASA_FFT_FIB_LIG_CELL_RESULTS_FINAL_4R.csv")
PATH.data.FFT.CN = file.path(PATH.data.FFT, "NASA_FFT_Project_CN_Data_4R.csv")
PATH.data.FFT.SLA_LMA = file.path(PATH.data.FFT, "NASA_FFT_SLA_LMA_Data_v2_4R_updated_new.csv")
PATH.results = "FFT_means.Rdata"
PATH.spectra = "~/Documents/Dropbox/FFT_spectra/FFT_spectra_unique.csv"
PATH.info = "~/Documents/"

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

#### Load results
load("run_results/FFTdata.Rdata")

## Basic statistics on individuals
leaf.means <- all.r[,lapply(.SD, mean), by = Spectra]
leaf.sd <- all.r[,lapply(.SD, sd), by=Spectra]

## Log statistics on individuals
all.r.log <- data.table(Spectra = all.r$Spectra)
for(cn in names(all.r)[2:5]) all.r.log[, cn := log(all.r[,cn, with=F]), with=F]
leaf.logmean <- all.r.log[,lapply(.SD, mean), by=Spectra]
leaf.logsd <- all.r.log[,lapply(.SD, mean), by=Spectra]

rename <- function(dt, s)
        setnames(dt, names(dt)[-1], sprintf("%s.%s", names(dt)[-1], s))

rename(leaf.means, "mu")
rename(leaf.sd, "sd")
rename(leaf.logmean, "lmu")
rename(leaf.logsd, "lsd")

## Merge data tables
setkey(leaf.means, Spectra)
setkey(leaf.sd, Spectra)
setkey(leaf.logmean, Spectra)
setkey(leaf.logsd, Spectra)
fft.data <- leaf.means[leaf.sd][leaf.logmean][leaf.logsd]
info.orig <- fread(PATH.spectra, header=TRUE, select = c(1:9,12,13))
setkey(info.orig, Spectra)
results <- fft.data[info.orig]

## Average duplicates
FFT2 <- FFT.all[, lapply(.SD, mean, na.rm=TRUE), by=c("Sample_Name", "Sample_Year")]

## Merge results with traits
fftdat <- merge(x=results, y=FFT2, by=mergeby.lower, all=TRUE)

### Sort out missing values
sample.name.regex <- "^([A-Za-z]+)([0-9]+[B]*)_([A-Za-z]+)_([BMTS]+)([ANO23]{0,1})[SAMP2]*"
fftdat[is.na(Species), Site := gsub(sample.name.regex, "\\1", Sample_Name)]
fftdat[is.na(Species), Plot := gsub(sample.name.regex, "\\1\\2", Sample_Name)]
fftdat[is.na(Species), Age := gsub(sample.name.regex, "\\5", Sample_Name)]
fftdat[is.na(Species), Height := gsub(sample.name.regex, "\\4", Sample_Name)]
fftdat[is.na(Species), Species := gsub(sample.name.regex, "\\3", Sample_Name)]
fftdat[Age=="T", Age := NA]
fftdat[Age=="0", Age := "N"]
fftdat[Age=="A", Age := "O"]
fftdat[is.na(Age) | Age == "", Age := "broadleaf"]
PATH.speciesinfo <- file.path("~/Documents", "Dropbox",
                              "NASA_TE_PEcAn-RTM_Project", "Data",
                              "FFT_species_info_csv.csv")
species.info <- fread(PATH.speciesinfo, header=TRUE)
species.info[,Scientific := sprintf("%s %s", Genus, Species)]
species.info[, Species := NULL]
setnames(fftdat, "Species", "Label")
setkey(fftdat, "Label")
setkey(species.info, "Label")
uk <- unique(c(fftdat[,Label], species.info[,Label]))
fft.full <- fftdat[species.info[J(uk)]]
fft.spec <- fft.full[!is.na(N.mu)]

## Fix species missing information

#fftdat["BEGL", list("Scientific
# exclude <- c("ANGE", "TYLA") ## As per Shawn's recommendation
save(fft.full, fft.spec, file="data/FFT_full.Rdata")
