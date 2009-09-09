###############################################################################
# function: findStretch()
# description: finds the stretch in the logical vector where the percentage of ones 
#              >= to the given percentage.
###############################################################################
findStretch <- function(x, stretch=5, percentage=0.5)
{
    if(sum(x) == 0) #no significant's at all
        return(NULL)   
    
    stretch <- stretch - 1	
    stretches <- data.frame()
    nx <- length(x)	
    
    y <- which(x[1:(nx-stretch)] == 1) #significant start positions
    ny <- length(y)
    
    i <- 1
    while(i  <= ny){
        if(sum(x[y[i]:(y[i]+stretch)]) < (stretch + 1)*percentage){
            i <- i + 1
            next
        }	
        j <- 1
        while(y[i]+j < nx - stretch & x[y[i]+stretch+j]) 
            j <- j + 1
        #store significant stretch
        start <- y[i]; end <- (y[i]+stretch+j-1)
        score <- 100*sum(x[start:end])/(end-start+1)
        stretches <- rbind(stretches, c(start, end, score))
        
        if(y[i]+2*stretch+j > y[ny]) #we are not at the end and another stretch fits
            break
        i <- which(y > y[i]+stretch+j)[1]		
    }
    
    if(nrow(stretches) == 0)
        return(NULL)
    colnames(stretches) <- c("start", "end", "score")	
    stretches	
}

###############################################################################
# function: findConsecutive()
# description: finds the consecutive regions
###############################################################################
findConsecutive <- function(xlogical)
{	
    if(sum(xlogical) == 0) #no significant's at all
         return(NULL)   
    x <- seq_len(length(xlogical))
    idx <- which(xlogical)
    breaks <- c(0, which(diff(idx) != 1), length(idx))			
    nbreaks <- length(breaks)
    
    consecutive <- data.frame(start=x[idx[breaks[1:(nbreaks-1)] + 1]],
            end=x[idx[breaks[2:nbreaks]]],
            score=rep(100, nbreaks-1))	
    consecutive
}

###############################################################################
# function: findWindow()
# description: finds the windows in the logical vector where the percentage of ones 
#              >= to the given percentage.
###############################################################################
findWindow <- function(x, y, window=1e2, percentage=0.5)
{
    if(sum(y) == 0) #no significant's at all
        return(NULL)   
    
    z <- withinWindow(x=x, y=x, w=window)	
    score <- apply(z, 1, function(zi){ tmp <- y[do.call(":", as.list(zi))]; sum(tmp)/length(tmp)})
    window <- data.frame(start=z[,1], end=z[,2], score=score)
    window[score >= percentage, ]
}

###############################################################################
# function: plotCytobands()
# description: add the Cytobands of a chromosome to an existing plot, taken from
#              SNPchip.
###############################################################################
plotCytobands <- function(chromosome, xlim, ylim, new=TRUE, label.cytoband=TRUE,
        cex.axis=1, outer=FALSE, taper=0.15,	db="homo_sapiens_core_40_36b", ...)
{
    
    #cytoband <- sim.update.chrom.table(db)
    data("chrom.table", package="SIM")
    
    table <- chrom.table[chrom.table$chr == chromosome,]	
    rownames(table) <- paste(table$arm, table$band, sep="")	
    table.p <- table[table$arm == "p", ]
    table.q <- table[table$arm == "q", ]	
    p.bands <- nrow(table.p)
    
    if(missing(ylim))
        ylim <- c(0.5, 1)	
    if(missing(xlim))
        xlim <- c(table$start[1], table$end[nrow(table)])
    
    cut.left <- c()
    cut.right <- c()
    for (i in 1:nrow(table)) {
        if (i == 1) {
            cut.left[i] <- TRUE
            cut.right[i] <- FALSE
        }
        else if (i == p.bands) {
            cut.left[i] <- FALSE
            cut.right[i] <- TRUE
        }
        else if (i == (p.bands + 1)) {
            cut.left[i] <- TRUE
            cut.right[i] <- FALSE
        }
        else if (i == nrow(table)) {
            cut.left[i] <- FALSE
            cut.right[i] <- TRUE
        }
        else {
            cut.left[i] <- FALSE
            cut.right[i] <- FALSE
        }
    }
    
    for (i in 1:nrow(table)) {
        if (as.character(table[i, "stain"]) == "stalk") {
            cut.right[i - 1] <- TRUE
            cut.left[i] <- NA
            cut.right[i] <- NA
            cut.left[i + 1] <- TRUE
        }
    }
    
    include <- table[, "end"] > xlim[1] & table[, "start"] < xlim[2]
    
    table <- table[include, ]
    
    n <- nrow(table)
    table[1, "start"] <- max(xlim[1], table[1, "start"])
    table[n, "end"] <- min(xlim[2], table[n, "end"])
    
    
    cut.left <- cut.left[include]
    cut.right <- cut.right[include]
    if (new) {
        xx <- c(0, table[n, "end"])		
        plot(xx, ylim, xlim=xlim, type="n", xlab="", ylab="",	axes=TRUE, yaxs="i", ylim=c(ylim[1], ylim[2] + 1))
    }
    bot <- ylim[1]
    top <- ylim[2]	
    h <- top - bot
    p <- taper
    
    getStain <- function(stain) {
        switch(stain, 
                gneg="grey100", 
                gpos25="grey90",
                gpos50="grey70", 
                gpos75="grey40", 
                gpos100="grey0",
                gvar="grey100", 
                stalk="brown3", 
                acen="brown4",
                "white")
    }
    
    color <- sapply(table$stain, getStain)
    
    for (i in 1:n) {
        start <- table[i, "start"]
        last <- table[i, "end"]
        delta <- (last - start)/4
        
        
        if(is.na(cut.left[i]) & is.na(cut.right[i])) {
            delta <- (last - start)/3
            segments(start + delta, bot, start + delta, top)
            segments(last - delta, bot, last - delta, top)
        }
        else if (cut.left[i] & cut.right[i]) {
            yy <- c(bot + p * h, bot, bot, bot + p * h, top - p * h, top, top, top - p * h)
            polygon(c(start, start + delta, last - delta, last, last, last - delta, start + delta, start), yy, col=color[i])
        }
        else if (cut.left[i]) {
            yy <- c(bot + p * h, bot, bot, top, top, top - p * h)
            polygon(c(start, start + delta, last, last, start + delta, start), yy, col=color[i])
        }
        else if (cut.right[i]) {
            yy <- c(bot, bot, bot + p * h, top - p * h, top, top)
            polygon(c(start, last - delta, last, last, last - delta, start), yy, col=color[i])
        }
        else {
            polygon(c(start, last, last, start), c(bot, bot, top, top), col=color[i])
        }
    }
    
    my.x <- (table[, "start"] + table[, "end"])/2
    if(label.cytoband)
        axis(1, at=my.x, labels=rownames(table), las=2, line=NA)
    
}


###############################################################################
# function: sim.plot.overlapping.indep.dep.features()
# description: plots and summarize the dependent and independent data together.
# TODO: overlapping not yet shown
###############################################################################
sim.plot.overlapping.indep.dep.features <- function(input.regions, 
        input.region.indep=NULL,
        adjust.method="BY", 
        log=FALSE, 
        significance=0.2,		
        max.pow=5,
        z.threshold= c(-3,3),
        summarize=c("consecutive", "stretch", "window"),
        stretch=10,
        window=1e6,
        percentage=0.5,
        xlim=NULL,
        pdf=FALSE,
        method=c("full", "smooth", "window", "overlap"),	
        run.name="analysis_results", ...)
{
    method <- match.arg(method)
    summarize <- match.arg(summarize)	
    
    #check significance
    if(!is.numeric(significance) | (significance < 0) | (significance > 1))
        stop("'significance' parameter should be a numeric value between 0 and 1! The input is: ", significance)
    
    input.regions <- unlist(sapply(input.regions, predefinedRegions, USE.NAMES=FALSE))
    
    ofPDF <- file.path(run.name, 
            paste("OverlappingPlot", method, adjust.method, format(significance, scientific = TRUE, digits=1), ".pdf", sep=""))
    
    cat("Producing overlapping plot ...\n")
    
    table.dep <- tabulate.top.dep.features(input.regions=input.regions, adjust.method=adjust.method, method=method, 
            significance=1, run.name=run.name)
    
    table.indep <- tabulate.top.indep.features(input.regions=input.regions, input.region.indep=input.region.indep, method=method, 
            adjust.method=adjust.method, significance=significance, decreasing=TRUE, 
            z.threshold=c(0, 0), run.name=run.name)
    
    if(pdf){
        options(error=dev.off)
        pdf(ofPDF, ...)
        
        on.exit({ dev.off()
                    cat("... overlapping plot stored with file name: ", ofPDF, "\n", sep="")
                    options(error=NULL)})
        
    } else {
        opar <- par(no.readonly=TRUE)
        on.exit(par(opar))
    }
    
    for(input.region in input.regions)
    {	
        region.dep <- getGenomicRegion(input.region)
        
        cat("... input region:", input.region, "\n")
        
        tableDep <- table.dep[[input.region]]		
        tableIndep <- table.indep[[input.region]]
        
        rescale <- 1e9
        chr <- as.numeric(region.dep$chr)
        
        x.dep <- tableDep$absolute.start.position - rescale*chr		
        x.indep <- tableIndep$absolute.start.position - rescale*chr
        
        tableDep <- tableDep[order(x.dep),]
        tableIndep <- tableIndep[order(x.indep),]
        
        x.dep <- x.dep[order(x.dep)]
        x.indep <- x.indep[order(x.indep)]
        
        y.dep <- tableDep$P.values        
        y.indep <- tableIndep$mean.Z.scores
        
        layout(matrix(c(1,2,3)), heights=c(2,2,1))
        par(mar=c(1.1, 4.1, 4.1, 3.1))
        
        main <- paste("Results for region", input.region, "according to location")		 
        xlab <- "" 
        ylab <- "P-value"		
        ylim <- c(0, 1)
        
        #maybe some checks on xlim
        if(is.null(xlim))
            xlim <- range(c(x.dep, x.indep))
        
        x <- switch(summarize, 
                consecutive=findConsecutive(y.dep <= significance),
                stretch=findStretch(y.dep <= significance, stretch=stretch, percentage=percentage),
                window=findWindow(x.dep, y.dep <= significance, window=window, percentage=percentage))
        
        ylog <- ifelse(log, "y", "")
        ##if log == TRUE do data transformation	
        if(ylog == "y") {
            log.sign <- 1/significance	
            ylim <- c(1, 10^max.pow)
            y.dep <- 1/y.dep
            ylab <- expression(-log[10]("P-value"))	
        }
                
        plot.new()		
        plot.window(xlim=xlim, ylim=ylim, log=ylog, xaxs="i")
        
        border <- "#5F9EA080"
        col <- "#5F9EA080"
        
        if(!is.null(x)) {
            for(i in 1:nrow(x))
                rect(x.dep[x$start[i]], ylim[1], x.dep[x$end[i]], ylim[2], col=col, border=border)
        }
        
        points(x.dep, y.dep, col=1, pch=16, cex=0.5)
                
        if(ylog == "y") {			
            abline(h=log.sign, col=1, lty=2)
            axis.at <- 10^c(0:max.pow)
            ##draw the major tick marks and label them using plotmath
            axis(2, at=axis.at, tcl=-1,	labels=parse(text=paste("10^-", 0:max.pow, sep="")))			
            ##now do the minor ticks, at 1/10 of each power of 10 interval
            axis(2, at=1:10 * rep(axis.at[-1] / 10, each=10), tcl=-0.5, labels=FALSE)
            axis(4, at=log.sign, labels=significance, col.ticks=1, col.axis= 1, las=2)
        } else {			
            abline(h=significance, col=1, lty=2)
            axis(2)
            axis(4, at=significance, labels=significance, col.ticks=1, col.axis= 1, las=2)
        }
        
        box()
        title(main=main, xlab=xlab, ylab=ylab)	
        
        x.dep.sign.pos <- matrix(c(x.dep[x$start], x.dep[x$end]), ncol=2, byrow=TRUE)
        
        par(mar=c(5.1, 4.1, 1.0, 3.1))
        
        ##mean zscore plot					
        xlab <- "base pair position (Mb)" 
        ylab <- "mean Z-scores"
        ylim <- range(c(y.indep[is.finite(y.indep)], z.threshold))
        
        plot.new()	
        plot.window(xlim=xlim, ylim=ylim, xaxs="i")
        
        col <- "#CD333380"
        border <- "#CD333380"
        
        x <- switch(summarize, 
                consecutive=findConsecutive(y.indep <= z.threshold[1] | y.indep >= z.threshold[2]),
                stretch=findStretch(y.indep <= z.threshold[1] | y.indep >= z.threshold[2], stretch=stretch, percentage=percentage),
                window=findWindow(x.indep, y.indep <= z.threshold[1] | y.indep >= z.threshold[2], window=window, percentage=percentage))
        
        if(!is.null(x)) {
            for(i in 1:nrow(x))
                rect(x.indep[x$start[i]], ylim[1], x.indep[x$end[i]], ylim[2], col=col, border=border)
        }		
        
        
        points(x.indep, y.indep, col=1, pch=16, cex=0.5)		
        abline(h=z.threshold, col=1, lty=2)
        axis(2)		
        axis(4, at=z.threshold[1], labels=z.threshold[1], col.ticks=1, col.axis=1, las=2)
        axis(4, at=z.threshold[2], labels=z.threshold[2], col.ticks=1, col.axis=1, las=2)
        
        Mb <- 1000000
        ax <- axTicks(1)		
        axis(1, at=ax, labels=as.double(ax/Mb))
        box()
        title(main="", xlab=xlab, ylab=ylab)	
        
        x.indep.sign.pos <- matrix(c(x.indep[x$start], x.indep[x$end]), ncol=2, byrow=TRUE)
        
        ##plot 3
        #par(mar=c(5.1, 4.1, 1.0, 3.1))
        #ylim <- c(0,1)
        #plot.new() 
        #plot.window(xlim=xlim, ylim=ylim, xaxs="i")
        
        #get the overlap
               
        
        ##plot chrom.table for the chromosome		
        par(mar=c(5, 4.1, 0.1, 3.1))	
        plot.new()	
        plot.window(xlim=xlim, ylim=c(-1, 1), xaxs="i")
        plotCytobands(chr, xlim=xlim, ylim=c(-0.5, 0.5), new=FALSE)
        
        xlim <- NULL
    }	
}
