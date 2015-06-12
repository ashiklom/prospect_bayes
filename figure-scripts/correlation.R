#' Correlation analysis figure
#' Author: Alexey Shiklomanov

# {{{ Head
if(grepl("figure-scripts", getwd())) setwd("..")
source("preprocess.fft.R")

library(ggplot2)
library(gridExtra)
source("preprocess.fft.R")
source("figure-scripts/common.R")
#source("figure-scripts/ggplot.statsmooth.R")
# }}}

# {{{ Variable definitions
pvars <- c('N.mu', 'Cab.mu', 'Cw.mu', 'Cm.mu')
cor.vars <- c("Perc_N", "Perc_C", "CNRatio",
                  "EWT_g_cm2", "LMA_g_DW_cm2", "SAMPLE_dN15", "ADF_PERC_DW",
                  "ADL_PERC_DW", "CELL_PERC_DW")
fft.cor.vars <- c(pvars, cor.vars)
fft.pairs.vars <- c(fft.cor.vars, "succession", "is.early", 'is.late')

#' Labeller for parameters
pvars.label <- list( "N.mu" = "N", "Cab.mu" = "Cab", "Cw.mu" = "Cw", "Cm.mu" = "Cm")
cor.vars.label <- list('Perc_N' = 'Perc N', 'Perc_C' = 'Perc C', 'CNRatio' = 'CN Ratio',
    'EWT_g_cm2' = 'EWT', 'LMA_g_DW_cm2' = 'LMA', 'SAMPLE_dN15' = 'dN15',
    'ADF_PERC_DW' = 'ADF', 'ADL_PERC_DW' = 'ADL', 'CELL_PERC_DW' = 'Cellulose')
vars.labeller <- function(variable, value) {
    if(variable == "cor.var") return(cor.vars.label[value])
    else if(variable == "pvar") return(pvars.label[value])
}
# }}}

fft[, is.early := succession == "early"][, is.late := succession == "late"]
fft.h[, is.early := succession == "early"][, is.late := succession == "late"]
fft.c[, is.early := succession == "early"][, is.late := succession == "late"]

# {{{ Panel function
#' Determine which correlations are significant
panel.lmsig <- function(x, y, f1, f2, f3, ...){
    dat <- data.frame(x, y)
    dat$f1 <- as.factor(f1)
    dat$f2 <- as.factor(f2)
    dat$f3 <- as.factor(f3)
    modlist <- list()
    modlist[["0"]] <- lm(y ~ 1, data=dat)
    modlist[["1"]] <- lm(y ~ x, data=dat)
    for(f in c("f1", "f2", "f3")){
        for(s in c("+", "*")){
            fs <- paste(f, s, sep=" ")
            form <- sprintf("y ~ %s x", fs)
            modlist[[fs]] <- lm(form, data=dat)
        }
    }
    mod.aic <- sapply(modlist, AIC)
    best <- which(mod.aic == min(mod.aic))
    print(best)
    bestf <- gsub(" [+*]", '', names(best))
    bestmod <- modlist[[best]]
    bestsum <- summary(bestmod)$coefficients
    print(rownames(bestsum))
    bestslopes <- bestsum[grepl("x", rownames(bestsum)),1]
    print(bestslopes)
    if(length(bestslopes) < 1) bestslopes <- 0
    bestints <- bestsum[!(rownames(bestsum) %in% names(bestslopes)),1]
    print(bestints)
    col.pal <- c("red", "green", "blue")
    if(bestf %in% c("0","1")){
        plot(dat[,"x"], dat[,"y"])
        abline(coef=c(bestints, bestslopes))
        print("line")
        return()
    }
    plot(dat[,"x"], dat[,"y"], col=col.pal[as.factor(dat[,bestf])])
    abline(coef=c(bestints[1], bestslopes[1]), col=col.pal[1])
    if(length(bestslopes) > 1){
        for(i in 2:length(bestslopes)){
            abline(coef=c(bestints[1] + bestints[i], bestslopes[1] + bestslopes[i]), col=col.pal[i])
            print("line")
        }
    } else if(length(bestints) > 1){
        for(i in 2:length(bestints)){
            abline(coef=c(bestints[1] + bestints[i], bestslopes), col=col.pal[i])
            print("line")
        }
    }
    print("")
}
# }}}
fft.hsub <- as.data.frame(fft.h[, fft.pairs.vars, with=F])
pdf("test.pdf", height=30, width=10)
layout(matrix(1:36, 9, 4, byrow=FALSE))
#par(mfrow=c(length(cor.vars), length(pvars)))
for(x in pvars){
    for(y in cor.vars){
        panel.lmsig(fft.hsub[,x], fft.hsub[,y], fft.hsub[, "succession"], fft.hsub[, "is.early"], fft.hsub[, "is.late"])
    }
}
dev.off()


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

sink("out.txt")
print(summary(mod1))
print(summary(mod3))
sink()
