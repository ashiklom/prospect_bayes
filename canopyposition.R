source("preprocess.fft.R")
library(ggplot2)
library(gridExtra)

## Canopy height
chp <- ggplot(fft.h) + aes(x=Height, fill=Height) + geom_violin() +
        xlab("Canopy position") +
        stat_summary(fun.y=mean, geom="point", size=2) + 
        stat_summary(fun.y=mean, geom="line", aes(group=1), size=1) +
        scale_x_discrete(labels=c("Bottom", "Middle", "Top")) +
        scale_fill_manual(values=c("lightgreen", "green", "darkgreen")) +
        guides(fill=FALSE) + 
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(size=11),
              axis.text = element_text(size=8))
ch.N <- chp + aes(y=N.mu) + ylab("Structure") + ggtitle("Hardwood")
ch.Cab <- chp + aes(y=Cab.mu) + ylab("Chlorophyll")
ch.Cw <- chp + aes(y=Cw.mu) + ylab("Water")
ch.Cm <- chp + aes(y=Cm.mu) + ylab("LMA")
ccp <- chp %+% fft.c + theme(axis.title.y = element_blank())
cc.N <- ccp + aes(y=N.mu) + ggtitle("Conifer")
cc.Cab <- ccp + aes(y=Cab.mu) 
cc.Cw <- ccp + aes(y=Cw.mu)
cc.Cm <- ccp + aes(y=Cm.mu)
png("manuscript/figures/canopy.png", width=6,height=6,units="in",res=300)
grid.arrange(ch.N, cc.N,
             ch.Cab, cc.Cab,
             ch.Cw, cc.Cw,
             ch.Cm, cc.Cm, ncol=2)
dev.off()
