tabulate.pvals <- function(input.regions="all chrs", 
        adjust.method="BY", 
        bins=c(0.001,0.005,0.01,0.025,0.05,0.075,0.10,0.20,1.0), 
        significance.idx=8, 
        order.by, 
        decreasing=TRUE, 
        method=c("full", "smooth", "window", "overlap"), 
        run.name="analysis_results")
{
    
    cat("Tabulate P-values ...\n")
    
    ##CHECK INPUT ARGUMENTS	
    method <- match.arg(method)
    
    if(!is.numeric(significance.idx) | significance.idx <= 0 | significance.idx > length(bins))
        stop("Wrong column selected for 'significance.idx': ", significance.idx, "\n")
    
    if(!missing(order.by))
        if(!(order.by %in% c("input.region", bins, "%")))
            stop("Wrong column given to order on, given column: ", order.by, "\n")
    
    ifPvaluesInputRegion <- file.path(run.name, method, "gpvals.pat.c.%s") 
    
    ##CHECK FILE EXISTENCE
    checkFiles(c(run.name, ifPvaluesInputRegion))	
    
    input.regions <- unlist(sapply(input.regions, predefinedRegions, USE.NAMES=FALSE))
    
    bins <- sort(bins)
    retval <- data.frame()
    for(input.region in input.regions)
    {	
        cat("... input region:", input.region, "\n")
        
        # Get the p-values for the current chromosome and adjust them for mulitple
        # testing.
        raw.pvals <- dget(sprintf(ifPvaluesInputRegion, input.region))		
        
        p.values <- p.adjust(raw.pvals, method=adjust.method)
        
        # Create the bins that. Each bin holds the number of p-values more or 
        # equally significant than its value (as specified by "bins").		
        comp.bins <- sapply(bins, function(x) sum(p.values <= x))
        
        # The last column holds the percentage of significant p-values on the chromosome.
        comp.bins <- c(comp.bins, round(100*(comp.bins[significance.idx] / length(p.values)), 2))
        
        # Cast comp.bins to an array and add the current set of bins as a row to the
        # matrix that we'll return.
        comp.bins <- matrix(comp.bins, nrow=1)
        retval <- rbind(retval, comp.bins)
    }
    
    retval <- cbind(input.regions, retval)
    
    # Cast the matrix to a data frame and set the colnames.
    retval <- as.data.frame(retval)
    colnames(retval) <- c("input.region", bins, "%")
    
    # Sort if necessary.
    if(!missing(order.by)){
        indices <- order(retval[,order.by], decreasing=decreasing)		
        retval <- retval[indices, ]
    }
    
    colnames(retval)[ncol(retval)] <- paste(significance.idx, "th (%)", sep="")
    
    return(retval)	
}

