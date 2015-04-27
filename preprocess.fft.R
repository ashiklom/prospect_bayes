library(data.table)
load("data/FFT_full.Rdata")

## Remove unknowns from analysis
bad.label <- c("UNKN1")
bad.pft <- NA
fft <- fft.spec[!(Label %in% bad.label)][!(PFT %in% bad.pft)]
fft[, PFT := tolower(PFT)]

# Combine North and South pines into single PFT
fft[PFT %in% c("south pine", "north pine"), PFT := "early conifer"]

# Add 'nontree' designator to shrubs and grass
fft[PFT == "shrub", PFT := "nontree shrub"]
fft[PFT == "grass", PFT := "nontree grass"]

## Pull out succession and plant type
fft.sp <- data.frame(do.call(rbind,
                             (strsplit(as.character(fft$PFT), " "))))
fft[, c("succession","plant.type") := fft.sp]

## Create ordered factors
succ.order <- c("nontree", "early", "mid", "late")
fft[, succession := factor(succession, levels=succ.order)]
type.order <- c("grass", "shrub", "conifer", "hardwood")
fft[, plant.type := factor(plant.type, levels=type.order)]
pft.order <- unique(fft[order(as.integer(succession), 
                              as.integer(plant.type)), PFT])
fft$PFT <- fft[, factor(PFT, levels=pft.order)]
label.order <- unique(fft[order(as.integer(PFT)), Label])
fft[, factor(Label, levels=label.order)]

## Set type-specific data.tables
fft.t <- fft[succession != "nontree"]
fft.h <- fft[plant.type == "hardwood"]
fft.c <- fft[plant.type == "conifer"]
