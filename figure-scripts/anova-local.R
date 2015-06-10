#' Complex analysis of variance for FFT data
# {{{ Head
if(grepl("figure-scripts", getwd())) setwd("..")

library(ggplot2)
library(gridExtra)
library(MASS)

source("preprocess.fft.R")
source("figure-scripts/common.R")
# }}}

# {{{ Support functions
reg.vars <- c('plant.type', 'succession', 'Label', 'Height',
              'Age', 'isnew', 'Site', 'Plot', 'Family')
p.vars <- c('N.mu', 'Cab.mu', 'Cw.mu', 'Cm.mu')

tofactor <- function(x){
    if(typeof(x) %in% c("character", "integer")) return(as.factor(x)) 
    else return(x)
}
datfix <- function(dat){
    dat[, isnew := (Age == "N")]
    dat <- dat[, c(p.vars, reg.vars), with=F]
    dat <- dat[, lapply(.SD, tofactor)]
    return(dat)
}
# }}}
# {{{ Fix data
fft <- fft[plant.type %in% c("hardwood", "conifer"),]
fft <- datfix(fft)
fft.c <- datfix(fft.c)
fft.h <- datfix(fft.h)
# }}}
# {{{ Formula setup
form.list <- c("a" = "%s ~ 1 + plant.type + succession + plant.type:succession + Label:plant.type:succession + Height + isnew + Site + Site:Plot",
               "h" = "%s ~ 1 + succession + Label:succession + Height + Site + Site:Plot",
               "c" = "%s ~ 1 + succession + Label:succession + Height + isnew + Site + Site:Plot")
form <- function(d,x) formula(sprintf(form.list[d], x))
# }}}

# {{{ Run models
mod.fft <- lapply(p.vars, function(x) lm(form('a',x), data=fft))
mod.fft.h <- lapply(p.vars, function(x) lm(form('h',x), data=fft.h))
mod.fft.c <- lapply(p.vars, function(x) lm(form('c',x), data=fft.c))
# }}}

# {{{ Dictionary and labeller
var.levels <- list('plant.type' = 'Plant type',
                   'succession' = 'Succession',
                   'plant.type:succession' = 'PFT',
                   'plant.type:succession:Label' = 'Species',
                   'succession:Label' = 'Species',
                   'Site' = 'Site',
                   'Site:Plot' = 'Plot',
                   'Height' = 'Canopy position',
                   'isnew' = 'Needle age',
                   'Residuals' = 'Residuals')
var.colors <- c('Plant type' = 'darkgreen',
                'Succession' = 'red',
                'PFT' = 'purple',
                'Species' = 'yellow',
                'Site' = 'brown',
                'Plot' = 'black',
                'Canopy position' = 'skyblue',
                'Needle age' = 'orange',
                'Residuals' = 'grey')
# }}}
# {{{ Themes
th.global <- theme_bw() +
    theme(text = element_text(size=11),
          axis.text = element_text(size=7),
          axis.title = element_text(size=9),
          legend.text = element_text(size=8),
          legend.title = element_text(size=10))
th.nox <- theme(axis.title.x = element_blank())
th.noy <- theme(axis.title.y = element_blank())
# }}}

# {{{ Plot function
anova.plot <- function(mod){
    require(ggplot2)
    require(reshape2)
    av.list <- lapply(mod, anova)
    vars.list <- rownames(av.list[[1]])
    ss.dat <- data.frame(sapply(av.list, function(x) x$'Sum Sq'/sum(x$'Sum Sq')))
    ss.dat$vars <- vars.list
    colnames(ss.dat) <- c(p.vars, "vars")
    ss.dat$vars <- sapply(ss.dat$vars, function(x) var.levels[[x]]) 
    cols <- var.colors
    ss.dat$vars <- factor(ss.dat$vars, levels=unlist(var.levels))
    colnames(ss.dat) <- c("N", "Cab", "Cw", "Cm", "vars")
    dat <- melt(ss.dat, id.vars="vars")
    plt <- ggplot(dat, aes(x=variable, y=value, fill=vars)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values = cols) +
        labs(x = "PROSPECT parameter", y = "Fraction of\nvariance explained")
    return(plt)
}
# }}}

# {{{ Generate plots
aov.all.l <- anova.plot(mod.fft) + th.global + 
    guides(fill = guide_legend(title="Predictor",
                               title.position="top",
                               title.hjust=0.5,
                               direction="horizontal",
                               nrow=2))
aov.all <- aov.all.l + th.global + guides(fill=FALSE) + th.nox + th.noy
aov.h <- anova.plot(mod.fft.h) + th.global + guides(fill=FALSE) + th.nox
aov.c <- anova.plot(mod.fft.c) + th.global + guides(fill=FALSE) + th.noy
all.legend <- get_legend(aov.all.l)
png.plot("manuscript/figures/variance-decomposition.png", h=4, w=4)
grid.arrange(arrangeGrob(aov.all, aov.h, aov.c, all.legend,
                         heights=c(1,1,1,0.5),ncol=1))
dev.off()
# }}}
