library(data.table)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(PEcAnRTM)
library(GGally)

## Set paths
path.db <- file.path("~", "Documents", "Dropbox")
path.fft <- file.path(path.db, "FFT_spectra", "FFT_spectra_unique.csv")

## Read FFT spectra
fft.spec <- fread(path.fft, header=TRUE, drop=c(10,11,14:21,22:71))

meltspec <- function(ss){
    mt <- melt(ss, measure.vars = sprintf("Wave_%d", 400:2500))
    mt$Wavelength <- as.numeric(gsub("Wave_([0-9]+)","\\1",mt$variable))
    mt$variable <- NULL
    return(mt)
}


## Observed spectra example
species.list <- c("ACRU", "QURU")
fft.melt <- meltspec(fft.spec[Species %in% species.list])
plt <- ggplot(fft.melt) + 
    aes(x=Wavelength, y=value, color=Species, group=Sample_Name) +
    geom_line(size=0.17, alpha=0.6) +
    labs(x = "Wavelength (nm)", y = "Reflectance") +
    guides(color = guide_legend(override.aes=list(alpha=1, size=2)))
plot(plt)

# Observed spectra -- png
png("Figures/observed.png", height=6, width=8, units="in", res=300)
plot(plt)
dev.off()


## Sample inversion
obs <- fft.spec[Species == "ACRU", 12:2112, with=FALSE]
obs <- t(as.matrix(obs))[,1]
ss <- invert_prospect(obs, 20000, c(1.8, 25, 0.01, 0.01))

bt <- seq(1000, 20000, by=15)
samps <- ss[bt,]

par(mfrow=c(2,2))
for(i in 1:4) plot(samps[,i], type='l', xlab="", ylab="")
sdf <- data.frame(samps[,1:4])

pairs.theme <- theme_bw() + 
        theme(text = element_text(size=30),
              axis.text = element_text(size=rel(0.4)),
              axis.title = element_text(size=rel(0.75)),
              plot.margin = unit(c(1,1,1,1), "lines"))
theme_set(pairs.theme)
blank <- grid.rect(gp=gpar(col="white"))
pair <- ggplot(sdf) + geom_point() 
dens <- ggplot(sdf) + geom_density() + 
        theme(axis.title = element_blank())
N.Cab <- pair + aes(x=N, y=Cab) + ylab("Chlorophyll") + 
        theme(axis.title.x = element_blank())
N.Cw <- pair + aes(x=N, y=Cw) + ylab("Water") + 
        theme(axis.title.x = element_blank())
N.Cm <- pair + aes(x=N, y=Cm) + ylab("LMA") + xlab("Structure")
Cab.Cw <- pair + aes(x=Cab, y=Cw) + 
        theme(axis.title = element_blank())
Cab.Cm <- pair + aes(x=Cab, y=Cm) + xlab("Chlorophyll") + 
        theme(axis.title.y = element_blank())
Cw.Cm <- pair + aes(x=Cw, y=Cm) + xlab("Water") + 
        theme(axis.title.y = element_blank())
Nd <- dens + aes(x=N)
Cabd <- dens + aes(x=Cab)
Cwd <- dens + aes(x=Cw)
Cmd <- dens + aes(x=Cm) + coord_flip()
png("Figures/pairs.png", width=14.2, height=8.5, units="in", res=300)
grid.arrange(Nd, blank, blank, blank,
             N.Cab, Cabd, blank, blank,
             N.Cw, Cab.Cw, Cwd, blank,
             N.Cm, Cab.Cm, Cw.Cm, Cmd, nrow=4, ncol=4)
dev.off()
                                     
rs <- 1:200
wl <- 400:2500
cp <- colorRampPalette(c("red", "blue"))
colors <- cp(length(rs))
colors <- adjustcolor(colors, alpha=0.6)
png("Figures/inversion.png", height=3, width=8, units="in", res=300)
par(mfrow=c(1,2), mai=c(0.45, 0.45, 0.1, 0.1))
plot(wl, rowMeans(obs), type='l', col="black", lwd=3, ylim=c(0,0.5),
     xlab='', ylab='')
plot(wl, prospect(4, ss[rs[1],]), type='l', col=colors[1], ylim=c(0,0.5),
     xlab='', ylab='')
for(r in 2:length(rs)) 
    lines(wl, prospect(4, ss[rs[r],]), col=colors[r])
dev.off()


## Setup parameter plots 
source("preprocess.fft.R")
succ.colors <- scale_color_manual(values = c("green3", "blue", "red"))

## Water content
water.base <- ggplot() + aes(x=Cw.mu, y=EWT_g_cm2) + 
        geom_point(size=3) + 
        geom_abline(linetype="dashed") + 
        xlim(0,0.12) + ylim(0,0.04) +
        theme(axis.title = element_blank(),
              legend.title = element_text(16),
              legend.text = element_text(12),
              legend.key.height = unit(1.5, "lines"),
              plot.margin = unit(c(1,1,1,1), "lines"))
water.all <- water.base %+% fft.t + aes(color=plant.type) +
        ggtitle("All trees") + 
        theme(legend.position=c(0.65, 0.13)) +
        guides(color = guide_legend(title="Plant type"))
water.h <- water.base %+% fft.h + aes(color=succession) +
        ggtitle("Hardwood") + 
        theme(legend.position=c(0.65, 0.13)) +
        xlim(0,0.025) + ylim(0,0.02) +
        guides(color = guide_legend(title="Succession")) +
        succ.colors
water.c <- water.h %+% fft.c + aes(color=succession) +
        xlim(0,0.12) + ylim(0,0.04) +
        ggtitle("Conifer")
png("Figures/water.png", height=7.76, width=13.63, units="in", res=300)
grid.arrange(water.all, water.h, water.c,  ncol=3)
dev.off()

## LMA
lma.base <- ggplot() + aes(x=Cm.mu, y=LMA_g_DW_cm2) + 
        geom_point(size=3) + 
        geom_abline(linetype="dashed") + 
        theme(axis.title = element_blank(),
              legend.title = element_text(size=28),
              legend.text = element_text(size=24),
              legend.key.height = unit(1.4, "lines"),
              plot.margin = unit(c(1,1,1,1), "lines"))
lma.all <- lma.base %+% fft.t + aes(color=plant.type) +
        ggtitle("All trees") +
        theme(legend.position=c(0.65, 0.13)) +
        guides(color = guide_legend(title="Plant type"))
lma.h <- lma.base %+% fft.h + aes(color=succession) +
        ggtitle("Hardwood") + 
        theme(legend.position=c(0.68, 0.12)) +
        guides(color = guide_legend(title="Succession")) +
        succ.colors
lma.c <- lma.h %+% fft.c + 
        ggtitle("Conifer") + 
        theme(legend.position=c(0.68, 0.12),
              axis.text.y = element_blank())
png("Figures/lma.png", height=7.76, width=13.63, units="in", res=300)
grid.arrange(lma.all, lma.h, lma.c, ncol=3)
dev.off()

lma.hin <- ggplot(fft.h) + 
        aes(x=Cm.mu, y=LMA_g_DW_cm2, color=succession) +
        geom_point(size=1) +
        geom_abline(linetype="dashed") +
        xlim(0,0.015) + ylim(0,0.015) +
        theme(text = element_text(size = 12),
              axis.title = element_blank(),
              plot.margin = unit(rep(0,4), "lines")) +
        guides(color = FALSE) + succ.colors
png("Figures/lma_inset.png",height=1.9,width=2.2,units="in", res=300)
plot(lma.hin)
dev.off()

## Canopy height
chp <- ggplot(fft.h) + aes(x=Height, fill=Height) + geom_violin() +
        xlab("Canopy position") +
        stat_summary(fun.y=mean, geom="point", size=5) + 
        stat_summary(fun.y=mean, geom="line", aes(group=1), size=2) +
        scale_x_discrete(labels=c("Bottom", "Middle", "Top")) +
        scale_fill_manual(values=c("lightgreen", "green", "darkgreen")) +
        guides(fill=FALSE) + 
        theme(axis.title.x = element_blank(),
              axis.text = element_text(size=18))
ch.N <- chp + aes(y=N.mu) + ylab("Structure") + ggtitle("Hardwood")
ch.Cab <- chp + aes(y=Cab.mu) + ylab("Chlorophyll")
ch.Cw <- chp + aes(y=Cw.mu) + ylab("Water")
ch.Cm <- chp + aes(y=Cm.mu) + ylab("LMA")
ccp <- chp %+% fft.c + theme(axis.title.y = element_blank())
cc.N <- ccp + aes(y=N.mu) + ggtitle("Conifer")
cc.Cab <- ccp + aes(y=Cab.mu) 
cc.Cw <- ccp + aes(y=Cw.mu)
cc.Cm <- ccp + aes(y=Cm.mu)
png("Figures/canopyheight.png", width=15.5,height=14,units="in",res=300)
grid.arrange(ch.N, cc.N,
             ch.Cab, cc.Cab,
             ch.Cw, cc.Cw,
             ch.Cm, cc.Cm, ncol=2)
dev.off()
