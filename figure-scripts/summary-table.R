#' Summary table

library(xtable)
if(grepl("figure-scripts", getwd())) setwd("..")
source("preprocess.fft.R")
source("figure-scripts/common.R")

vars <- c('N', 'Cab', 'Cw', 'Cm')
vars.mu <- paste0(vars, '.mu')

# {{{ Summarize function
summarize <- function(datn, succ=TRUE){
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
    m <- formatC(as.matrix(qmed[,vars.mu, with=FALSE]), format="f", digits=2)
    ql <- formatC(as.matrix(qlow[,vars.mu, with=FALSE]), format="f", digits=2)
    qh <- formatC(as.matrix(qhigh[,vars.mu,with=FALSE]), format="f", digits=2)
    sm <- matrix(sprintf("%s (%s, %s)", m, ql, qh), nrow=nrow(m), ncol=ncol(m))
    sm.dat <- data.table(sm)
    setnames(sm.dat, vars)
    if(succ) {
        sm.dat[, succession := qlow[,succession]]
        setkey(sm.dat, "succession")
    }
    sm.dat[, dat := datn]
    return(sm.dat)
}
# }}}

fft.h[, Cw.mu := 1000 * Cw.mu][,Cm.mu := 10000 * Cm.mu]
fft.c[, Cw.mu := 1000 * Cw.mu][,Cm.mu := 10000 * Cm.mu]

summary.h <- summarize("fft.h", T)
summary.c <- summarize("fft.c", T)

sum.hc <- rbind(summary.h, summary.c)
sum.hc[dat=="fft.c", dat := "Conifer"]
sum.hc[dat=="fft.h", dat := "Hardwood"]
varnames <- c("N",
              "Cab (\\si{\\micro\\gram\\per\\square\\centi\\meter})",
              "Cw (\\si{\\centi\\metre} $\\times 10^3$)",
              "Cm (\\si{\\gram\\per\\square\\meter})")
sum.hc[c(1,3,4,6), dat := ""]
setcolorder(sum.hc, c("dat","succession", vars))
setnames(sum.hc, c("dat", "succession", vars), c("Plant type", "Succession", varnames))

out.caption <- "
Summary of inversion estiamtes of PROSPECT parameters.
Values are reported as mean (95\\% CI)
"
out.caption <- gsub("\\n", " ", out.caption)
out.tab <- xtable(sum.hc, caption=out.caption, label="tab:summary")

out.tab.pre <- print(out.tab, file="", 
      include.rownames=F,
      sanitize.colnames.function=identity, 
      size="small")
out.tab.post <- gsub("(\\{\\\\small)", "\\\\centerline\\{\\1", out.tab.pre)
out.tab.post <- gsub("(\\\\end\\{tabular\\}\\n)", "\\1\\}", out.tab.post)

cat(out.tab.post, file="manuscript/tables/summary-table.tex")
