library(data.table)
load("data/FFT_full.RData")

## Fix needle age
fft.full[Age=="T", Age := NA]
fft.full[Age=="0", Age := "N"]
fft.full[Age=="A", Age := "O"]
fft.full[is.na(Age) | Age == "", Age := "broadleaf"]

## Remove unknowns from analysis
bad.label <- c("UNKN1")
bad.pft <- NA
fft.f <- fft.full[!(Label %in% bad.label)][!(PFT %in% bad.pft)]
fft.f[, PFT := tolower(PFT)]

# Combine North and South pines into single PFT
fft.f[PFT %in% c("south pine", "north pine"), PFT := "early conifer"]

# Add 'nontree' designator to shrubs and grass
fft.f[PFT == "shrub", PFT := "nontree shrub"]
fft.f[PFT == "grass", PFT := "nontree grass"]

## Pull out succession and plant type
fft.sp <- data.frame(do.call(rbind,
                             (strsplit(as.character(fft.f$PFT), " "))))
fft.f[, c("succession","plant.type") := fft.sp]

## Create ordered factors
succ.order <- c("nontree", "early", "mid", "late")
fft.f[, succession := factor(succession, levels=succ.order)]
type.order <- c("grass", "shrub", "conifer", "hardwood")
fft.f[, plant.type := factor(plant.type, levels=type.order)]
pft.order <- unique(fft.f[order(as.integer(succession), 
                              as.integer(plant.type)), PFT])
fft.f$PFT <- fft.f[, factor(PFT, levels=pft.order)]
label.order <- unique(fft.f[order(as.integer(PFT)), Label])
fft.f[, factor(Label, levels=label.order)]

## Set type-specific data.tables
fft.f <- fft.f[succession != "nontree"]
fft.h <- fft.f[plant.type == "hardwood"]
fft.c <- fft.f[plant.type == "conifer"]

## Save
save(fft.f, fft.h, fft.c, file="data/FFT.processed.RData")
