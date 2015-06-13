#' Correlation figure
# {{{ Source and library
if(grepl("figure-scripts", getwd())) setwd("..")
source("preprocess.fft.R")
source("figure-scripts/common.R")

library(ggplot2)
library(reshape2)
library(gridExtra)
# }}}

# {{{ Variable definitions
pvars <- c('N.mu', 'Cab.mu', 'Cw.mu', 'Cm.mu')
cor.vars <- c("Perc_N", "Perc_C", "CNRatio",
                  "EWT_g_cm2", "LMA_g_DW_cm2", "SAMPLE_dN15", "ADF_PERC_DW",
                  "ADL_PERC_DW", "CELL_PERC_DW")
succ.levels <- c("early", "mid", "late")
p.label <- list( "N.mu" = "N", "Cab.mu" = "Cab", "Cw.mu" = "Cw", "Cm.mu" = "Cm")
v.label <- list('Perc_N' = 'Perc N', 'Perc_C' = 'Perc C', 'CNRatio' = 'CN Ratio',
    'EWT_g_cm2' = 'EWT', 'LMA_g_DW_cm2' = 'LMA', 'SAMPLE_dN15' = 'dN15',
    'ADF_PERC_DW' = 'ADF', 'ADL_PERC_DW' = 'ADL', 'CELL_PERC_DW' = 'Cellulose')
col.pal <- c("green2", "blue", "red")
# }}}

# {{{ Plot function
plot.lm <- function(dat, pp){
    par(pp)
    pp <- par(mfcol=c(length(cor.vars), length(pvars)), 
          mar=c(0.2,0.2,0.2,0.2), cex.axis=0.8)
    dat[, succession := factor(succession, levels=succ.levels)]
    for(p in pvars){
        for(v in cor.vars){
            col1 <- p == pvars[1]
            rown <- v == cor.vars[length(cor.vars)]
            pv <- paste(p, v)
            f.basic <- sprintf("%s ~ %s", v, p)
            f.succ <- sprintf("%s ~ succession * %s", v, p)
            mod.basic <- lm(f.basic, data=dat)
            mod.succ <- lm(f.succ, data=dat)
            aic.basic <- AIC(mod.basic)
            aic.succ <- AIC(mod.succ)
            if(aic.succ < aic.basic){
                modsum <- summary(mod.succ)
                if(col1 && rown)
                    plot(x=dat[[p]], y=dat[[v]], col=col.pal[dat[,succession]],
                         pch=20, cex=0.2)
                else if (col1)
                    plot(x=dat[[p]], y=dat[[v]], col=col.pal[dat[,succession]],
                         pch=20, cex=0.2, xaxt='n')
                else if (rown)
                    plot(x=dat[[p]], y=dat[[v]], col=col.pal[dat[,succession]],
                         pch=20, cex=0.2, yaxt='n')
                else
                    plot(x=dat[[p]], y=dat[[v]], col=col.pal[dat[,succession]],
                         pch=20, cex=0.2, xaxt='n', yaxt='n')
                for(s in succ.levels){
                    dat.sub <- dat[succession==s,]
                    mod.sub <- lm(f.basic, data=dat.sub)
                    cf <- summary(mod.sub)$coefficients
                    if(all(cf[,4] < 0.05)) abline(coef=cf[,1], col=col.pal[match(s,succ.levels)])
                }
            } else {
                modsum <- summary(mod.basic)
                if(col1 && rown)
                    plot(x=dat[[p]], y=dat[[v]], 
                         pch=20, cex=0.2)
                else if (col1)
                    plot(x=dat[[p]], y=dat[[v]], 
                         pch=20, cex=0.2, xaxt='n')
                else if (rown)
                    plot(x=dat[[p]], y=dat[[v]], 
                         pch=20, cex=0.2, yaxt='n')
                else
                    plot(x=dat[[p]], y=dat[[v]], 
                         pch=20, cex=0.2, xaxt='n', yaxt='n')
                cf <- modsum$coefficients
                if(all(cf[,4] < 0.05)) abline(coef=cf[,1])
            }
            if(col1) mtext(v.label[[v]], side=2, cex=0.7, line=2)
            if(rown) mtext(p.label[[p]], side=1, cex=0.7, line=2)
        }
    }
}
# }}}
# {{{ Draw plots
pp0 <- list(omi=c(0.4, 0.4, 0.7, 0.1), cex.main=2)
leg <- function() {
    pp <<- par(usr=c(0,1,0,1), xpd=NA)
    legend(-2.1, 10.5, succ.levels, 
           col=col.pal, lty=1, horiz=TRUE)
}
png.plot("manuscript/figures/correlation-hardwood.png", h=6, w=4)
pp <- pp0
plot.lm(fft.h, pp)
title("Hardwood", outer=TRUE, line=3)
leg()
dev.off()
png.plot("manuscript/figures/correlation-conifer.png", h=6, w=4)
pp <- pp0
plot.lm(fft.c[EWT_g_cm2 < 0.1,], pp)
title("Conifer", outer=TRUE, line=3)
leg()
dev.off()
# }}}

