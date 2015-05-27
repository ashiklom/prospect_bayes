source("preprocess.fft.R")

## Plot
library(ggplot2)
library(gridExtra)

succ.colors <- scale_color_manual(values = c("green3", "blue", "red"))

## Water content
water.base <- ggplot() + aes(x=Cw.mu, y=EWT_g_cm2) + 
        geom_point(size=1.5) + 
        geom_abline(linetype="dashed") + 
        theme(axis.title = element_blank(),
              axis.text = element_text(size=8),
              legend.title = element_text(size=11),
              legend.text = element_text(size=9),
              legend.key.height = unit(0.7, "lines"),
              legend.key.width = unit(0.7, "lines"),
              plot.margin = unit(c(0.4,0.4,0.4,0.4), "lines"))
water.all <- water.base %+% fft.t + aes(color=plant.type) +
        theme(legend.position="bottom") +
        scale_x_continuous(breaks=c(0, 0.03, 0.06, 0.09, 0.12)) +
        xlim(0,0.125) + ylim(0,0.04) +
        guides(color = guide_legend(title="Plant type",
                                    title.position="top"))
water.h <- water.base %+% fft.h + aes(color=succession) +
        theme(legend.position="bottom") +
        xlim(0,0.028) + ylim(0,0.02) +
        scale_x_continuous(breaks=c(0, 0.005, 0.01, 0.015, 0.02, 0.025)) +
        guides(color = guide_legend(title="Succession",
                                    title.position="top")) +
        succ.colors
water.c <- water.h %+% fft.c + aes(color=succession) +
        xlim(0,0.125) + ylim(0,0.04) +
        scale_x_continuous(breaks=c(0, 0.03, 0.06, 0.09, 0.12)) +
        guides(color = FALSE) + 
        theme(plot.margin = unit(c(0.4,0.4,3.2,0.4), "lines"))
png("manuscript/figures/water.png", height=4, width=6, units="in", res=300)
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
png("manuscript/figures/lma.png", height=7.76, width=13.63, units="in", res=300)
grid.arrange(lma.all, lma.h, lma.c, ncol=3)
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
