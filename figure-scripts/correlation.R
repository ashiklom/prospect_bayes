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
# }}}

# {{{ Function to extract values of interest
get.lm <- function(dat){
    dalist <- list()
    dslist <- list()
    for(p in pvars){
        for(v in cor.vars){
            pv <- paste(p, v)
            f.basic <- sprintf("%s ~ %s", v, p)
            f.succ <- sprintf("%s ~ succession * %s", v, p)
            dalist[[pv]] <- lm(f.basic, data=dat)
            dslist[[pv]] <- lm(f.succ, data=dat)
        }
    }
    dasum <- lapply(dalist, summary)
    dssum <- lapply(dslist, summary)
    dacoefs <- data.frame(lapply(dasum, '[[', 'coefficients'))
    dscoefs <- data.frame(lapply(dssum, '[[', 'coefficients'))
    dacoefs$varname <- rownames(dacoefs)
    dscoefs$varname <- rownames(dscoefs)
    dacoefs$varname <- gsub("N\\.mu", "x", dacoefs$varname)
    dscoefs$varname <- gsub("N\\.mu", "x", dscoefs$varname)
    dscoefs$varname <- gsub("(.*)", "succ_\\1", dscoefs$varname)
    damelt <- melt(dacoefs, id.vars="varname")
    dsmelt <- melt(dscoefs, id.vars="varname") 
    var.rxp <- "((N|Cab|Cw|Cm)\\.mu)\\.(.*?)\\.+(.*)"
    damelt$pvar <- gsub(var.rxp, "\\2", damelt$variable)
    damelt$cvar <- gsub(var.rxp, "\\3", damelt$variable)
    damelt$pc <- gsub(var.rxp, "\\1 \\3", damelt$variable)
    damelt$param <- gsub(var.rxp, "\\4", damelt$variable)
    dsmelt$pvar <- gsub(var.rxp, "\\2", dsmelt$variable)
    dsmelt$cvar <- gsub(var.rxp, "\\3", dsmelt$variable)
    dsmelt$pc <- gsub(var.rxp, "\\1 \\3", dsmelt$variable)
    dsmelt$param <- gsub(var.rxp, "\\4", dsmelt$variable)
    dacast <- data.table(dcast(damelt, varname + pvar + cvar + pc ~ param))
    dscast <- data.table(dcast(dsmelt, varname + pvar + cvar + pc ~ param))

    dar2 <- sapply(dasum, '[[', 'r.squared')
    dar2.dt <- data.table(rsquare=dar2, index=names(dar2))
    daaic <- sapply(dalist, AIC)
    daaic.dt <- data.table(aic=daaic, index=names(daaic))
    setkey(dar2.dt, index)
    setkey(daaic.dt, index)
    setkey(dacast, pc)
    dacast <- dacast[dar2.dt][daaic.dt]

    dsr2 <- sapply(dssum, '[[', 'r.squared')
    dsr2.dt <- data.table(rsquare=dsr2, index=names(dsr2))
    dsaic <- sapply(dslist, AIC)
    dsaic.dt <- data.table(aic=dsaic, index=names(dsaic))
    setkey(dsr2.dt, index)
    setkey(dsaic.dt, index)
    setkey(dscast, pc)
    dscast <- dscast[dsr2.dt][dsaic.dt]

    dat <- rbind(dacast, dscast)
    dat[grep("succ", varname), modtype := "succession"]
    dat[is.na(modtype), modtype := "basic"]

    setnames(dat, c("Pr...t..", "Std..Error"), c("p.value", "SE"))

    return(dat)
}
# }}}

hdat <- get.lm(fft.h)
hdat[, plant.type := "hardwood"]
cdat <- get.lm(fft.c)
cdat[, plant.type := "conifer"]
dat <- rbind(hdat, cdat)
datm <- melt(dat1, measure.vars=c('Estimate', 'p.value', 'SE', 't.value', 'rsquare', 'aic'))
datc <- dcast(datm, pvar + cvar + pc + plant.type ~ varname + variable)

alpha <- 0.05
dat[p.value < 0.05, rsquare2 := rsquare]

plt <- ggplot(dat, aes(x=pvar, y=cvar, fill=rsquare2)) + 
    geom_tile(color="black") + 
    facet_grid(plant.type ~ modtype) +
    scale_fill_gradient(low="white", high="red")
plot(plt)







#' Labeller for parameters
pvars.label <- list( "N.mu" = "N", "Cab.mu" = "Cab", "Cw.mu" = "Cw", "Cm.mu" = "Cm")
cor.vars.label <- list('Perc_N' = 'Perc N', 'Perc_C' = 'Perc C', 'CNRatio' = 'CN Ratio',
    'EWT_g_cm2' = 'EWT', 'LMA_g_DW_cm2' = 'LMA', 'SAMPLE_dN15' = 'dN15',
    'ADF_PERC_DW' = 'ADF', 'ADL_PERC_DW' = 'ADL', 'CELL_PERC_DW' = 'Cellulose')
vars.labeller <- function(variable, value) {
    if(variable == "cor.var") return(cor.vars.label[value])
    else if(variable == "pvar") return(pvars.label[value])
}


fft.cor.vars <- c(pvars, cor.vars)
cor.plot <- function(dat){
    fft.cor <- dat[, fft.cor.vars, with=F]
    cor.all <- data.table(cor(fft.cor, use="pairwise.complete.obs"))
    cor.all$var1 <- colnames(cor.all)
    setcolorder(cor.all, c("var1", fft.cor.vars))
    cor.melt <- melt(cor.all, id.vars="var1")
    cor.melt <- cor.melt[variable %in% pvars,][!(var1 %in% pvars),]
    c1 <- ggplot(cor.melt) + aes(x=variable, y=var1, fill=value) +
        geom_tile() + 
        scale_fill_gradient2(low="red", mid="white", high="green") +
            xlab("Prospect parameter") +
            ylab("Measured variable")
    return(c1)
}

cor.all <- cor.plot(fft) + ggtitle("All")
cor.h <- cor.plot(fft.h) + ggtitle("Hardwood")
cor.c <- cor.plot(fft.c) + ggtitle("Conifer")

pairs.plot <- function(dat){
    dat.m1 <- melt(dat, measure.vars=pvars, variable.name="pvar", value.name="pval")
    dat.m2 <- melt(dat.m1, measure.vars=cor.vars, variable.name="cor.var", value.name="cor.val")
    plt <- ggplot(dat.m2) + 
        aes(x=pval, y=cor.val) +
        geom_point(size=0.8) +
        facet_grid(cor.var ~ pvar, scales="free", labeller=vars.labeller) +
        stat_smooth(method="lm", se=FALSE, fullrange=FALSE, size=1)
    return(plt)
}
#
theme.pairs <- theme_bw() + 
    theme(text = element_text(size=10),
          axis.text = element_text(size=rel(0.6)),
          strip.text = element_text(size=rel(0.7)),
          axis.title.y = element_blank())
fft.c.sub <- fft.c[EWT_g_cm2 < 0.1,][N.mu < 5,]
pairs.h <- pairs.plot(fft.h) + aes(color=succession) + xlab("Hardwood") + theme.pairs + succ.colors +
    guides(color=FALSE)
pairs.c <- pairs.plot(fft.c.sub) + aes(color=succession) + xlab("Conifer") + theme.pairs + succ.colors
png.plot(fname="manuscript/figures/correlation.png", h=7, w=9)
grid.arrange(arrangeGrob(pairs.h, pairs.c, nrow=1, widths=c(1,1.3)))
dev.off()

