#' Summary plot for all inversion results

library(ggplot2)
library(gridExtra)
setwd("..")
source("preprocess.fft.R")
source("figure-scripts/common.R")

#' Summary plot -- global theme
th.summary <- theme_bw() + 
    theme(text = element_text(size=10),
          axis.title = element_text(size=rel(1)),
          axis.text = element_text(size=rel(0.6))
          )
ww <- c(1, 1, 1)
#
pN <- ggplot() + aes(x=PFT, y=N.mu) + geom_violin() +
    ylab("Mesophyll structure") +
    th.summary +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
pN.h <- pN %+% fft.h + aes(fill=succession) + 
    succ.colors + 
    theme(axis.title.y = element_text(size=rel(1))) +
    guides(fill=FALSE)
pN.c <- pN %+% fft.c + aes(fill=succession) + 
    succ.colors +
    guides(fill=FALSE) +
    theme(axis.title = element_blank())
pN.nt <- pN %+% fft[succession=="nontree",] + geom_violin(fill="grey") +
    theme(axis.title = element_blank())
grobN <- arrangeGrob(pN.h, pN.c, pN.nt, nrow=1, widths=ww)
#
pCab.h <- pN.h + aes(y=Cab.mu) + ylab("Chlorophyll (ug/cm2)")
pCab.c <- pN.c + aes(y=Cab.mu)
pCab.nt <- pN.nt + aes(y=Cab.mu)
grobCab <- arrangeGrob(pCab.h, pCab.c, pCab.nt, nrow=1, widths=ww)
#
pCw.h <- pCab.h + aes(y=Cw.mu) + ylab("EWT (g/cm2)")
pCw.c <- pCab.c + aes(y=Cw.mu)
pCw.nt <- pCab.nt + aes(y=Cw.mu)
grobCw <- arrangeGrob(pCw.h, pCw.c, pCw.nt, nrow=1, widths=ww)
#
no.xlabs <- theme(axis.text.x = element_blank())
no.y <- theme(axis.title.y = element_blank())
xlab.nt <- scale_x_discrete("Non-tree", labels=c("Grass", "Shrub"))
pCm.h <- pCab.h + aes(y=Cm.mu) + ylab("LMA (g/cm2)") + xlab("Hardwoods") + th.summary + no.xlabs
pCm.c <- pCab.c + aes(y=Cm.mu) + xlab("Conifers") + th.summary + no.xlabs + no.y
pCm.nt <- pCab.nt + aes(y=Cm.mu) + th.summary + no.y + xlab.nt
grobCm <- arrangeGrob(pCm.h, pCm.c, pCm.nt, nrow=1, widths=ww)
#
leg.gg <- pN.h + guides(fill="legend") +
    theme(legend.direction="horizontal",
          legend.text=element_text(size=rel(0.8)),
          legend.title=element_text(size=rel(1))
          )
leg <- ggplot_gtable(ggplot_build(leg.gg))$grobs[[8]]
png.plot("manuscript/figures/summary.png", h=8, w=6)
grobs <- arrangeGrob(grobN, grobCab, grobCw, grobCm, leg, nrow=5, ncol=1,
                     heights=c(1,1,1,1.1,0.2))
grid.arrange(grobs)
dev.off()
