#' Summary table

library(xtable)
if(grepl("figure-scripts", getwd())) setwd("..")
source("preprocess.fft.R")
source("figure-scripts/common.R")

vars <- c('N', 'Cab', 'Cw', 'Cm')
vars.mu <- paste0(vars, '.mu')

summarize <- function(datn, succ=TRUE){# [--
    dat <- get(datn)
    if(succ){
        qlow <- dat[, lapply(.SD, quantile, 0.025, na.rm=TRUE), by=succession, .SDcols=vars.mu]
        qmed <- dat[, lapply(.SD, quantile, 0.5, na.rm=TRUE), by=succession, .SDcols=vars.mu]
        qhigh <- dat[, lapply(.SD, quantile, 0.975, na.rm=TRUE), by=succession, .SDcols=vars.mu]
    } else {
        qlow <- dat[, lapply(.SD, quantile, 0.025, na.rm=TRUE), .SDcols=vars.mu]
        qmed <- dat[, lapply(.SD, quantile, 0.5, na.rm=TRUE), .SDcols=vars.mu]
        qhigh <- dat[, lapply(.SD, quantile, 0.975, na.rm=TRUE), .SDcols=vars.mu]
    }
    m <- formatC(as.matrix(qmed[,vars.mu, with=FALSE]), format="f", digits=3)
    ql <- formatC(as.matrix(qlow[,vars.mu, with=FALSE]), format="f", digits=3)
    qh <- formatC(as.matrix(qhigh[,vars.mu,with=FALSE]), format="f", digits=3)
    sm <- matrix(sprintf("%s (%s, %s)", m, ql, qh), nrow=nrow(m), ncol=ncol(m))
    sm.dat <- data.table(sm)
    setnames(sm.dat, vars)
    if(succ) {
        sm.dat[, succession := qlow[,succession]]
        setkey(sm.dat, "succession")
    }
    sm.dat[, dat := datn]
    return(sm.dat)
}# --]

summary.h2 <- summarize("fft.h", F)
summary.c2 <- summarize("fft.c", F)

sum2 <- rbind(summary.h2, summary.c2)
sum2[dat=="fft.c", dat := "Conifer"]
sum2[dat=="fft.h", dat := "Hardwood"]
setcolorder(sum2, c("dat", vars))
setnames(sum2, "dat", "Plant type")


out.caption <- "
Summary of inversion estiamtes of PROSPECT parameters.
Values are reported as median (95\\% CI)
"
out.caption <- gsub("\\n", " ", out.caption)
out.tab <- xtable(sum2, caption=out.caption, label="tab:summary")

print(out.tab, file="manuscript/tables/summary-table.tex", include.rownames=F)
