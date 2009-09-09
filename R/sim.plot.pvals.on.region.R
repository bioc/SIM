###############################################################################
# function: sim.plot.pvals.on.region()
# description: make several diagnostic plots of the P-values and multiple testing
#              correction.
###############################################################################
sim.plot.pvals.on.region <-	function(input.regions=c("all chrs"),
        significance=0.2,	
        adjust.method="BY", 
        method=c("full", "smooth", "window", "overlap"), 
        run.name="analysis_results", ...)
{
    cat("Producing region plot ...\n")
    
    ##CHECK INPUT	
    method <- match.arg(method)
    
    if(!is.numeric(significance) | (any(significance < 0)) || (any(significance > 1)))		
        stop("'significance' parameter should be a numeric value between 0 and 1! The input is: ", paste(significance, collapse=", "))
    
    #dget	 
    ifDepPosInputRegion  <- file.path(run.name, method, "start.input.region.%s")
    ifPvalues <- file.path(run.name, method, "gpvals.pat.c.%s")
    
    ofPDF <- file.path(run.name, "pvalue.plots", paste("PvalueOnRegion", method, adjust.method, format(significance, scientific = TRUE, digits=1), ".pdf", sep=""))
    
    ##CHECK FILE EXISTENCE
    checkFiles(c(run.name, ifDepPosInputRegion, ifPvalues))
    
    # Open the graphics device.		
    pdf(ofPDF, ...)
    
    options(error=dev.off)
    
    input.regions <- unlist(sapply(input.regions, predefinedRegions, USE.NAMES=FALSE))
    
    for(input.region in input.regions)
    {
        cat("... input region:", input.region, "\n")	
        
        raw.pvals <- dget(sprintf(ifPvalues, input.region))		
        abs.start <- dget(sprintf(ifDepPosInputRegion, input.region))
        
        p.values <- p.adjust(raw.pvals, method=adjust.method)
        col <- "blue"
        pch <- 20
        cex <- ifelse(length(abs.start) < 1000, 1, 0.5)
        
        layout(rbind(c(1,2,3), c(4,4,4)), heights=c(2,2), respect=TRUE)
        # Histogram.
        par(mar=c(5.1, 4.1, 4.1, 0.1))
        hist(raw.pvals, main="raw P-values", xlab="P-value", col=col, xlim=c(0, 1))
        
        # Plot of raw p-values.		
        par(mar=c(5.1, 4.1, 4.1, 0.1))
        plot(sort(raw.pvals), ylim=c(0,1), main="raw P-values", xlab="genes", ylab="sorted p-values", pch=pch, col=col, cex=cex)		
        abline(a=0, b=1/length(raw.pvals), lty="dashed")
        
        # Plot of adjusted p-values.	
        par(mar=c(5.1, 4.1, 4.1, 0.1))
        plot(sort(p.values), ylim=c(0,1), main=paste(adjust.method, "-corrected\nP-values"), xlab="genes", 
                ylab="sorted adjusted p-values", pch=pch, col=col, cex=cex)			
        abline(h=significance, lty="dashed")
        
        # Plot the pvalues along the chromosome.
        par(mar=c(5.1, 4.1, 2.1, 0.1))
        plot(abs.start, p.values, main=paste("Input region:", input.region), ylim=c(0,1), cex=cex, pch=pch, col=col, 
                xaxt="n", xlab="absolute start position (Mb)", ylab=paste(adjust.method, "-corrected P-values"))
        Mb <- 1000000
        ax <- axTicks(1)		
        axis(1, at=ax, labels=as.integer(ax)/Mb)		
        abline(h= significance, lty="dashed")		
    }
    
    dev.off()
    options(error=NULL)
    cat("... region plot stored with file name: \n", ofPDF, "\n", sep="")
}

