###############################################################################
# function: sim.plot.pvals.on.genome()
# description: plots chromosome overview of the human genome and highlights
#              the significant P-values of the input region given by the user.
###############################################################################
sim.plot.pvals.on.genome <- function(input.regions="all chrs", 
        significance=c(0.05, 0.20), 
        adjust.method="BY", 
        method=c("full", "smooth", "window", "overlap"),		
        run.name="analysis_results",
        pdf=TRUE, 
        main="Significantly associated features",		
        ylab="Chromosomes",		
        ann=par("ann"),
        ...)
{
    data("chrom.table", package="SIM")
    
    ##CHECK INPUT	
    method <- match.arg(method)
    
    #check P-value significane parameter
    if(!is.numeric(significance) | (any(significance < 0)) || (any(significance > 1)))		
        stop("'significance' parameter should be a numeric value between 0 and 1! The input is: ", paste(significance, collapse=", "))
    
    #dget	 
    ifDepPosInputRegion <- file.path(run.name, method, "start.input.region.%s")
    ifPvalues <- file.path(run.name, method, "gpvals.pat.c.%s") 
    
    ofPDF <- file.path(run.name, "pvalue.plots",  
            paste("WholeGenomePlot", method, adjust.method, ".pdf", sep=""))
    
    ##CHECK FILE EXISTENCE
    checkFiles(c(run.name, ifDepPosInputRegion, ifPvalues))
    
    cat("Producing genome plot ...\n")
    
    pval <- significance[order(significance)]
    
    input.regions <- unlist(sapply(input.regions, predefinedRegions, USE.NAMES=FALSE))
    
    ##PREPARE PLOTTING PARAMETERS
    
    nchrs <- length(levels(chrom.table$chr))	
    
    chrsLength <- sapply(levels(chrom.table$chr), function(x) convertGenomicRegion(x)$end)
    
    xlines <- matrix(c(rep(0, nchrs), chrsLength), nrow=2, byrow=TRUE) 
    ylines <- matrix(c(nchrs:1, nchrs:1), nrow=2, byrow=TRUE)
    
    xpoints <- sapply(levels(chrom.table$chr), function(x) convertGenomicRegion(x, "p")$end)
    ypoints <- nchrs:1			
    
    hsegments <- 0.4
    
    #define colors
    col.chromosomes <- rep("black", nchrs)	
    col.centromers <- "purple"	
    col.legend <- c("#00ffff", "#0088ff", "#9f9f9f")
    
    ##PLOTTING		
    if(pdf){			
        options(error=dev.off)
        pdf(ofPDF, ...)
        
        on.exit({ dev.off()
                    cat("... overlapping plot stored with file name: \n", ofPDF, "\n", sep="")
                    options(error=NULL)})
        
    } else {
        opar <- par(no.readonly=TRUE)
        on.exit(par(opar))
    }
    
    par(mar=c(0.1, 4.1, 4.1, 0.1))
    
    plot.new()	
    plot.window(c(0, max(chrsLength)), c(0, nchrs))		
    matlines(xlines, ylines, col=col.chromosomes, lty=1)	
    
    for(input.region in input.regions)
    {
        region <- getGenomicRegion(input.region)
        
        cat("... input region:", input.region, "\n")
        
        chrom <- as.numeric(region$chr)
        
        raw.pvals <- dget(sprintf(ifPvalues, input.region))
        abs.start <- dget(sprintf(ifDepPosInputRegion, input.region))
        
        adjpval <- p.adjust(raw.pvals, method=adjust.method)
        
        xsegments <- matrix(c(abs.start, abs.start), nrow=2, byrow=TRUE)		
        ysegments <- matrix(rep(c(nchrs + 1 - chrom - hsegments, nchrs + 1 - chrom + hsegments), each=length(abs.start)), nrow=2, byrow=TRUE)
        
        indices <- which(adjpval > pval[2])
        if(length(indices) > 0)
            matlines(xsegments[,indices], ysegments[,indices], col=col.legend[3], lty=1, pch=1)
        indices <- which(adjpval < pval[2] & adjpval >= pval[1])
        if(length(indices) > 0)
            matlines(xsegments[,indices], ysegments[,indices], col=col.legend[2], lty=1, pch=1)
        indices <- which(adjpval <= pval[1])
        if(length(indices) > 0)
            matlines(xsegments[,indices], ysegments[,indices], col=col.legend[1], lty=1, pch=1)		
        lines(c(region$start, region$end), ylines[, chrom], col="orange", lty=1)
    }
    
    points(xpoints, ypoints, col=col.centromers, bg=col.centromers, pch=19)	
    axis(2, at=1:nchrs, labels=c("Y", "X", (nchrs-2):1), las=2) #new axis x, y i.s.o. 23, 24
    
    if(ann)
        title(main=main, ylab=ylab)	
    legendText1 <- substitute(p <= pMin, list(pMin=pval[1]))
    legendText2 <- substitute(paste(pMin < p, NULL <= pMax), list(pMin=pval[1], pMax=pval[2]))
    legendText3 <- substitute(p > pMax, list(pMax=pval[2]))	
    legend("right", do.call("expression", list(legendText1, legendText2, legendText3)), fill=col.legend, bty='n')
    
}

