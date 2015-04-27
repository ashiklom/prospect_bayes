if(grepl("Figures", getwd())) setwd("..")
source("preprocess.fft.R")
library(reshape2)

## ANOVA colors
colors.anova = scale_fill_manual(values = c("Canopy position"="red",
                                            "Species"="green",
                                            "Plant type"="darkgreen",
                                            "Succession"="orange",
                                            "Plot"="skyblue",
                                            "Site"="darkblue",
                                            "Residuals"="Grey"))

## ANOVAs
anova.var <- function(vd,vi, df){
    fm <- formula(sprintf("%s ~ %s", vd, paste(vi, collapse=" + ")))
    av <- anova(lm(fm, data=df))
    out <- av$"Sum Sq"/sum(av$"Sum Sq")
    names(out) <- rownames(av)
    return(out)
}
anova.df <- function(fft, f1){
    pv <- c("N.mu", "Cab.mu", "Cw.mu", "Cm.mu")
    av <- lapply(pv, anova.var, vi=f1, df=fft)
    dav <- data.frame(av)
    colnames(dav) <- pv
    dav$labels <- rownames(dav)
    dav <- data.table(dav)
    avm <- melt(dav, id.vars="labels")
    avm[labels == "Label", labels := "Species"]
    avm[labels == "Height", labels := "Canopy position"]
    avm[labels == "plant.type", labels := "Plant type"]
    avm[labels == "succession", labels := "Succession"]
    lvs <- c("Residuals", "Plot", "Site", "Species", 
             "Plant type", "Succession", "Canopy position")
    avm[, labels := factor(labels, lvs)]
    return(avm)
}
f1 <- c("plant.type", "succession", "Label", 
        "Height", "Site")
p.all <- anova.df(fft, f1)
p.tree <- anova.df(fft.t, f1)
p.h <- anova.df(fft.h, f1[-1])
p.c <- anova.df(fft.c, f1[-1])

av.plot <- ggplot(p.tree) + 
        aes(x=variable, y=value, fill=labels) +
        geom_bar(stat="identity") + 
        colors.anova +
        xlab("Variable") +
        ylab("Frac. variance explained") +
        labs(fill="\n\n\n\n\n\n\n\n\n\n\n\n") +
        scale_x_discrete(labels = c("Structure", "Chlorophyll",
                                    "Water", "LMA")) +
        theme(text = element_text(size=28),
              axis.title.y=element_text(size=24),
              axis.title.x=element_blank(),
              axis.text.y=element_text(size=16),
              axis.text.x=element_text(size=18,angle=-35, vjust=0.5),
              legend.text = element_text(size=26),
              legend.key.height = unit(2.5, "lines"),
              axis.title.x = element_blank())
av.all <- av.plot + guides(fill = guide_legend(ncol=3)) +
        ggtitle("All trees") + 
        theme(plot.margin=unit(c(1,1,0,1), "lines"))
png("Figures/anova.all.png", width=15.5, height=5.3, units="in", res=300)
plot(av.all)
dev.off()

av.hcp <- av.plot + guides(fill=FALSE) + 
        theme(plot.margin=unit(c(1,1,0,0), "lines"))
av.h <- av.hcp %+% p.h + ggtitle("Hardwood")
av.c <- av.hcp %+% p.c + ggtitle("Conifer") + 
        theme(axis.title.y = element_blank())
png("Figures/anova.pft.png", width=15.5,height=5.3,units="in",res=300)
grid.arrange(av.h, av.c, ncol=2)
dev.off()

## Mean with fixed effects
# require(rjags)
# jags.code <-" 
#     model{
#Priors
#         for(b in 1:nbeta){beta[b] ~ dnorm(0, 0.01)}
#         rinv ~ dgamma(0.01, 0.01)
# 
#Model
#         for(i in 1:nobs){
#             Ex[i] <- beta %*% X[i,]
#             ym[i] ~ dnorm(Ex[i], std[i])
#             y[i] ~ dnorm(Ex[i], rinv)
#         }
#     }
# "
# fft.mod <- fft[, list(Cw.mu, Cw.sd, plant.type, succession,
#                       Height, Age, Site, Plot)]
# fft.mod <- fft.mod[complete.cases(fft.mod)]
# obsmat <- model.matrix(Cw.mu ~ plant.type + succession + Height + 
#                        Age + Site + Plot, data=fft)
# jags.data <- list(y = fft$Cw.mu,
#                   std = fft$Cw.sd,
#                   X = obsmat,
#                   nbeta = ncol(obsmat),
#                   nobs = nrow(obsmat))
# model <- jags.model(file=textConnection(jags.code), 
#                     data=jags.data,
#                     n.chains=4)
# 
# samps <- coda.samples(model, n.iter=1000,
#                       variable.names = c("nbeta", "rinv"))
