library(data.table)
library(reshape2)

### Process SED files
i <- 1
read.SE <- function(path, count=TRUE){
        if(count){
                i <<- i+1
                print(i)
        }
        rawfile <- readLines(path)
        specstart <- grep("Data:", rawfile) + 2
        header <- rawfile[1:(specstart-4)]
        head.rxp <- "(.*?):(.*)"
        hnames <- gsub(head.rxp, "\\1", header)
        hvals <- gsub(head.rxp, "\\2", header)
        
        # Raw spectrum
        rawspec <- rawfile[specstart:length(rawfile)]
        speclist <- strsplit(rawspec, "\\t")
        speclist2 <- lapply(speclist, as.numeric)
        specdat <- do.call(rbind, speclist2)
        refl <- specdat[,4]
        wl <- specdat[,1]
        rlist <- c(as.list(hvals), as.list(refl))
        names(rlist) <- c(hnames, sprintf("Wave_%d", wl))
        return(rlist)
}

read.all.SE <- function(path, break.every = 500){
        se.rxp <- ".*\\.sed$"
        flist <- list.files(path, se.rxp, recursive = TRUE, full.names = TRUE)
        nf <- length(flist)
        nff <- length(flist) / break.every
        i <- 1
        biglist <- lapply(flist[1:break.every], read.SE)
        bigspec <- rbindlist(biglist)
        rm(biglist)
        for(n in nff){
                biglist <- lapply(flist[(n*break.every+1):min(((n+1)*break.every),nf)],
                                  read.SE)
                bigspec <- rbindlist(list(bigspec, biglist))
                rm(biglist)
        }        

        ## Process spec name
        
        
        return(bigspec)
}

