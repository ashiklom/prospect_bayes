#' Validation plot based on Water and LMA
source("../figure.common.R")

# Preprocess data {{{
# Keep only rows with inversion data
fft.full <- fft.f
fft.f <- fft.f[!is.na(N.mu)][sensor=="identity"]
fft.h <- fft.f[plant.type == "hardwood"]
fft.c <- fft.f[plant.type == "conifer"]

# Factor the sensor list
#fft.full[, sensor := factor(sensor, levels=sensor.list)]
# }}}

# Global theme {{{
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
# }}}

# Water content plot {{{
water.base <- ggplot() + aes(x=EWT_g_cm2, y=Cw.mu) + 
        geom_point(size=1) + 
        geom_abline(linetype="dashed") + 
        xlim(0,0.06) +
        xlab("Measured") +
        ylab("Inversion estimate") +
        th.global +
        theme(legend.position="right",
              legend.margin=unit(0, "lines"))
water.all <- water.base %+% fft.f + aes(color=plant.type) +
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
pdf("manuscript/figures/water.pdf", height=4, width=6)
grid.arrange(arrangeGrob(water.all),
             arrangeGrob(water.h, water.c, nrow=1, widths=c(1,1.4)),
             nrow=2)
dev.off()
# }}}

# LMA plot {{{
lma.base <- ggplot() + aes(x=LMA_g_DW_cm2, y=Cm.mu) + 
        geom_point(size=1) + 
        geom_abline(linetype="dashed") + 
        xlim(0,0.035) +
        xlab("Measured") +
        ylab("Inversion estimate") +
        th.global +
        theme(legend.position="right",
              legend.margin=unit(0, "lines"))
lma.all <- lma.base %+% fft.f + aes(color=plant.type) +
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
pdf("manuscript/figures/lma.pdf", height=4, width=6) 
grid.arrange(arrangeGrob(lma.all),
             arrangeGrob(lma.h, lma.c, nrow=1, widths=c(1,1.4)),
             ncol=1)
dev.off()
# }}}

# Error statistic table {{{
rmse <- function(mod, obs){
  nx <- min(sum(!is.na(mod)), sum(!is.na(obs)))
  rmse <- sqrt(sum((mod - obs)^2, na.rm=TRUE) / nx)
  bias <- mean(mod - obs, na.rm=TRUE)
  sepc <- sqrt(sum((mod - obs - bias)^2, na.rm=TRUE) / nx)
  cv <- sepc / mean(mod, na.rm=TRUE)
  rmspe <- sqrt(sum((mod/obs - 1)^2, na.rm=TRUE) / nx)
  out <- list(RMSE=rmse, BIAS=bias, SEPC=sepc, CV=cv, RMSPE=rmspe)
  return(out)
}
rmse.water <- fft.spec[, rmse(Cw.mu, EWT_g_cm2), by=plant.type][,param := "EWT"]
rmse.lma <- fft.spec[, rmse(Cm.mu, LMA_g_DW_cm2), by=plant.type][,param := "LMA"]
rmse.table <- rbind(rmse.water, rmse.lma)

setcolorder(rmse.table, c("param", "plant.type", "RMSE", "BIAS", "SEPC", "CV", "RMSPE"))
setnames(rmse.table, c("param", "plant.type"), c("Parameter", "Plant type"))
# }}}

# Sensor response plots {{{
#rmse.plot <- ggplot(rmse.melt) + aes(x=sensor, y=value, fill=plant.type) +
    #facet_wrap(stat~param, scales="free") + geom_bar(stat="identity", position="dodge") +
    #theme(text = element_text(size=6), axis.text.x = element_text(angle=90, hjust=1))
#png.plot("fft-error-by-sensor.png", 5, 5)
#plot(rmse.plot)
#dev.off()
# }}}

# Prepare xtable {{{
library(xtable)
cap <- "
Error statistics for equivalent water thickness (EWT) and leaf mass per unit area (LMA)
model estimates compared to inversion.
"
cap <- gsub("\\n", " ", cap)
out.tab <- xtable(rmse.table, caption=cap, digits=4, label="tab:water-lma")
out.tab.pre <- print(out.tab, file="", include.rownames=FALSE)
out.tab.post <- out.tab.pre
cat(out.tab.post, file="manuscript/tables/water-lma.tex")
# }}}

