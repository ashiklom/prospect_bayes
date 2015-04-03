## Exploratory figures
library(data.table)
library(ggplot2)
library(GGally)
load("data/FFT_full.Rdata")

## Remove unknowns from analysis
bad.label <- "UNKN1"
bad.pft <- c("Grass", "Shrub", NA)
fft <- fft.spec[Label != bad.label][PFT != bad.pft]


## Traits by PFT
pft.plot <- ggplot(fft) + aes(x=Height) + facet_wrap("PFT", scales = "free_y")
pft.plot + aes(y=N.mu) + stat_summary(fun.data="mean_cl_boot", geom="pointrange", aes(group=1)) +
        stat_summary(fun.data=median, geom="line", aes(group=1))
ggplot(fft.spec) + aes(x=PFT, y=N.mu, col=PFT) + geom_boxplot()
ggplot(fft.spec) + aes(x=PFT, y=Cab.mu, col=PFT) + geom_boxplot()
ggplot(fft.spec) + aes(x=PFT, y=Cw.mu, col=PFT) + geom_boxplot()
ggplot(fft.spec) + aes(x=PFT, y=Cm.mu, col=PFT) + geom_boxplot()

## Trends by canopy position
ggplot(fft.spec) + aes(x=PFT, y=N, col=Height) + geom_violin() + geom_line(stat="median")
ggplot(fft.spec) + aes(x=PFT, y=Cab, col=Height) + geom_violin()
ggplot(fft.spec) + aes(x=PFT, y=Cw, col=Height) + geom_violin()
ggplot(fft.spec) + aes(x=PFT, y=Cm, col=Height) + geom_violin()


## LMA vs. N
ggplot(fft.spec) + aes(x=N, y=LMA_g_DW_cm2, col=PFT) +
        geom_point() +
        xlim(1,3)

## EWT vs. Cw
ggplot(fft.spec) + aes(x=Cw, y=EWT_g_cm2, col=PFT) + 
        geom_point() +
        ylim(0,0.05)

## CN vs Cab
ggplot(fft.spec) + aes(x=Cab, y=CNRatio, col=PFT) + geom_point()

## Pairs plot of PROSPECT parameters
ggpairs(fft.spec, 12:15, color="PFT")
ggpairs(fft.spec[PFT != "North pine"], 12:15, color="PFT")
ggpairs(fft.spec[PFT != "North pine" & N < 4 & Cab < 120 & Cw <0.05 & Cm <0.05], 12:15, color="PFT")

