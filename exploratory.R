## Exploratory figures
library(data.table)
library(ggplot2)
library(GGally)
library(reshape2)
if(grepl("Figures", getwd())) setwd("..")
source("preprocess.fft.R")

# Setup plot themes and colors for consistency
theme_set(theme_bw(base_size=20))
colors.cbl = scale_color_manual(values = c("BL"="green4",
                                           "Con"="brown"))
colors.height = scale_color_manual(values = c("B"="red",
                                              "M"="yellow",
                                              "T"="blue"))
colors.pft = scale_color_manual(values = c("Early Hardwood"="red",
                                               "Mid Hardwood"="green",
                                               "Late hardwood"="blue"))
colors.anova = scale_fill_manual(values = c("Canopy position"="red",
                                             "Species"="green",
                                            "PFT"="darkgreen",
                                             "Plot"="skyblue",
                                             "Site"="darkblue",
                                             "Residuals"="Grey"))

## Water content
water <- ggplot(fft) + aes(x=Cw.mu, y=EWT_g_cm2) +
        geom_abline(intercept=0, slope=1, linetype="dashed") +
        xlab("Estimated from spectra") +
        ylab("Measured directly") +
        ggtitle("Equivalent water thickness (cm)")

water.all <- water + 
        ylim(0,0.04) + 
        aes(size=Cw.sd/Cw.mu, col=Con_vs_BL) + 
        geom_point()

water.byspec <- water + ylim(0,0.04) + facet_wrap("Label")

water.fit <- lm(EWT_g_cm2 ~ Cw.mu, data=fft.bl)
wb <- water.fit$coefficients[1]
wm <- water.fit$coefficients[2]

water %+% fft.bl + aes(col=Height) +
        geom_point() + 
        geom_abline(intercept=wb, slope=wm, col="black")


## Parameter covariance
ggpairs(fft, 5:8, color="Con_vs_BL", params=c(shape="+"))
ggpairs(fft.bl, 5:8, color="PFT")


## LMA vs. Cm
lma.cm <- ggplot(fft.bl) + aes(x=Cm.mu, y=LMA_g_DW_cm2, color=Con_vs_BL) + 
        geom_point() +
        xlab("Estimated from spectra") +
        ylab("Measured directly") +
        ggtitle("LMA (g/cm2)")
pft.lv <- c("Early hardwood", "Mid Hardwood", "Late hardwood")
lma.cm + aes(color=factor(PFT, levels=pft.lv)) + 
        geom_smooth(method="lm", se=F, aes(group=Label), fullrange=T) + 
        labs(color="") +
        geom_abline(linetype="dashed", size=2) + ylim(0,0.02)

## Canopy height
chp <- ggplot(fft.bl) + aes(x=Height) + geom_violin() +
        xlab("Canopy position") +
        stat_summary(fun.y=mean, geom="point") + 
        stat_summary(fun.y=mean, geom="line", aes(group=1))
chp + aes(y=N.mu) + ylab("Structural complexity")
chp + aes(y=Cab.mu) + ylab("Chlorophyll (ug/cm2)")
chp + aes(y=Cw.mu) + ylab("Water (cm)")
chp + aes(y=Cm.mu) + ylab("LMA (g/cm2)")


## Basic ANOVAs
av.N <- anova(lm(N.mu ~ Height + PFT + Label + Site + Plot, data=fft.bl))
av.Cab <- anova(lm(Cab.mu ~ Height + PFT + Label + PFT + Site + Plot, fft.bl))
av.Cw <- anova(lm(Cw.mu ~ Height + PFT + Label + Site + Plot, fft.bl))
av.Cm <- anova(lm(Cm.mu ~ Height + PFT + Label + Site + Plot, fft.bl))
av <- data.table(N = av.N$"Sum Sq",
                 Cab = av.Cab$"Sum Sq",
                 Cw = av.Cw$"Sum Sq",
                 Cm = av.Cm$"Sum Sq")
av2 <- av[,lapply(.SD, function(x) x/sum(x))]
av2$labels <- rownames(av.N)
avm <- melt(av2, id.vars="labels")
avm[labels == "Label", labels := "Species"]
avm[labels == "Height", labels := "Canopy position"]
lvs <- c("Residuals", "Plot", "Site", "Species", "PFT", "Canopy position")
ggplot(avm) + 
        aes(x=variable, y=value, fill=factor(labels, levels=lvs)) +
        geom_bar(stat="identity") + 
        colors.anova +
        xlab("Variable") +
        ylab("Fraction of variance explained") +
        labs(fill="")




