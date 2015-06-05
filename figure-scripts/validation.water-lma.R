#' Validation plot based on Water and LMA

setwd("..")
source("preprocess.fft.R")

## Plot
library(ggplot2)
library(gridExtra)

succ.colors <- scale_color_manual(values = c("green3", "blue", "red"))

## Global settings
th.global <- theme_bw() + 
    theme(text = element_text(size=10),
          axis.title = element_text(size=rel(1)),
          axis.text = element_text(size=rel(0.6)),
          legend.title = element_text(size=rel(1)),
          legend.text = element_text(size=rel(0.7)),
          legend.key.height = unit(0.7, "lines"),
          legend.key.width = unit(0.7, "lines"),
          plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines")
          )
png.plot <- function(fname, h=4, w=4, ...){
    png(fname, height=h, width=w, res=300, units="in", ...)
}

## Water content
water.base <- ggplot() + aes(x=EWT_g_cm2, y=Cw.mu) + 
        geom_point(size=1) + 
        geom_abline(linetype="dashed") + 
        xlim(0,0.06) +
        xlab("Measured") +
        ylab("Inversion estimate") +
        th.global +
        theme(legend.position="right",
              legend.margin=unit(0, "lines"))
water.all <- water.base %+% fft.t + aes(color=plant.type) +
        guides(color = guide_legend(title="Plant type",
                                    title.position="top")) +
        annotate("text", x=0, y=0.11, label="(a)", size=4)
water.h <- water.base %+% fft.h + aes(color=succession) +
        xlim(0,0.02) +
        guides(color=FALSE) +
        succ.colors +
        annotate("text", x=0.001, y=0.021, label="(b)", size=4)
water.c <- water.base %+% fft.c + aes(color=succession) +
        xlim(0,0.06) + succ.colors +
        annotate("text", x=0.003, y=0.11, label="(c)", size=4) +
        guides(color = guide_legend(title="Succession",
                                    title.position="top")) +
        theme(axis.title.y = element_blank())
png.plot("manuscript/figures/water.png", h=4, w=6)
grid.arrange(arrangeGrob(water.all),
             arrangeGrob(water.h, water.c, nrow=1, widths=c(1,1.4)),
             nrow=2)
dev.off()


## LMA
lma.base <- ggplot() + aes(x=LMA_g_DW_cm2, y=Cm.mu) + 
        geom_point(size=1) + 
        geom_abline(linetype="dashed") + 
        xlim(0,0.035) +
        xlab("Measured") +
        ylab("Inversion estimate") +
        th.global +
        theme(legend.position="right",
              legend.margin=unit(0, "lines"))
lma.all <- lma.base %+% fft.t + aes(color=plant.type) +
        guides(color = guide_legend(title="Plant type",
                                    title.position="top")) +
        annotate("text", x=0, y=0.082, label="(a)", size=4)
        #annotate("text", x=-Inf, y=Inf, label="(a)", size=4, hjust=1, vjust=1)
lma.h <- lma.base %+% fft.h + aes(color=succession) +
        xlim(0,0.0125) +
        guides(color=FALSE) +
        succ.colors +
        annotate("text", x=0.0005, y=0.015, label="(b)", size=4)
lma.c <- lma.base %+% fft.c + aes(color=succession) +
        xlim(0,0.035) +
        succ.colors +
        guides(color = guide_legend(title="Succession",
                                    title.position="top")) +
        annotate("text", x=0.002, y=0.08, label="(c)", size=4)
        theme(axis.title.y = element_blank())
png.plot("manuscript/figures/lma.png", h=4, w=6) 
grid.arrange(arrangeGrob(lma.all),
             arrangeGrob(lma.h, lma.c, nrow=1, widths=c(1,1.4)),
             ncol=1)
dev.off()

