# Sensor covariance plot
library(ggplot2)
library(PEcAnRTM)
library(gridExtra)
library(grid)
library(reshape2)

# Data identifier -- Each ID corresponds to a single set of true values
datnum <- 100    
parnames <- c("N", "Cab", "Car", "Cw", "Cm")
datlist <- list()
for(s in sensor.list){
    f.path <- sprintf("../data/some-simulations/%s.%d.RData", s, datnum)
    load(f.path)
    f.name <- gsub(".*/(.*)[.]RData", "\\1", f.path)
    f.list <- get(f.name)
    rm(list=f.name)
    sensor <- list(sensor = f.list$sensor)
    samples <- f.list$samples
    true.param <- as.list(f.list$true.param)
    names(true.param) <- sprintf("%s.true", parnames[1:5])
# Burn in and thin
    ngibbs <- nrow(samples)
    burnin <- floor(ngibbs/2)
    thin <- (ngibbs - burnin) / 5000
    bt <- seq(burnin, ngibbs, by=thin)
    samples.sub <- samples[bt,]
    samples.pairs <- samples.sub[,1:5]
    colnames(samples.pairs) <- parnames[1:5]
    samples.dat <- data.table(samples.pairs)[,sensor := s]
    datlist[[f.name]] <- samples.dat
}


# Create pairs plot
params <- c("N", "Cab", "Car", "Cw", "Cm")
names(true.param) <- params
param.pairs <- combn(params, 2, simplify=FALSE)
plot.sensor <- c("identity", "aviris.ng", "modis", "landsat8")
sensor.color <- c("orange", "blue", "dark green", "dark red")
names(sensor.color) <- plot.sensor
plot.data <- sampdat[sensor %in% s][, sensor := factor(sensor, levels=rev(plot.sensor))]
plot.all <- ggplot(plot.data) + aes(color=sensor, order=sensor) +
    geom_density2d(size=0.2) +
    stat_density2d(aes(fill=sensor), alpha=0.2,  geom="polygon") +
    scale_color_manual(values = sensor.color) + 
    scale_fill_manual(values = sensor.color) +
    guides(color = FALSE, fill = FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=5),
          axis.title = element_blank(),
          plot.margin = unit(c(0,0,0.05,0.05), "lines"))
plotlist <- list()
for(p in param.pairs){
    pc <- sprintf("%s-%s", p[1], p[2])
    py <- true.param[[p[2]]]
    px <- true.param[[p[1]]]
    plt <- plot.all + aes_string(x=p[1], y=p[2]) +
        geom_hline(y=py, linetype="dashed", color="black", size=0.6) +
        geom_vline(x=px, linetype="dashed", color="black", size=0.6)
    if(p[1] != "N") plt <- plt + theme(axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank())
    if(p[2] != "Cm") plt <- plt + theme(axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank())
    plotlist[[pc]] <- plt
}
blank <- grid.rect(gp = gpar(col="white"), draw=FALSE)
pdf(file="manuscript/figures/pairs-4.pdf", width=7, height=5)
grid.arrange(textGrob("N"), blank, blank, blank, blank,
             plotlist[["N-Cab"]], textGrob("Cab"), blank, blank, blank,
             plotlist[["N-Car"]], plotlist[["Cab-Car"]], textGrob("Car"), blank, blank,
             plotlist[["N-Cw"]], plotlist[["Cab-Cw"]], plotlist[["Car-Cw"]], textGrob("Cw"), blank,
             plotlist[["N-Cm"]], plotlist[["Cab-Cm"]], plotlist[["Car-Cm"]], plotlist[["Cw-Cm"]], textGrob("Cm"),
             nrow=5)
dev.off()

