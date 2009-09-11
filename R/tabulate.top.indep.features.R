tabulate.top.indep.features <- function(input.regions="all chrs", 
        input.region.indep=NULL,
        method=c("full", "smooth", "window", "overlap"), 
        adjust.method="BY", 
        significance=1, 
        decreasing=TRUE, 
        z.threshold=c(0, 0),				 
        run.name="analysis_results")
{ 
    cat("Tabulate results on independent data ...\n")
    
    ##CHECK INPUT ARGUMENTS	
    method <- match.arg(method)	
    
    #check significance
    if(!is.numeric(significance) | (significance < 0) | (significance > 1))
        stop("'significance' parameter should be a numeric value between 0 and 1! The input is: ", significance)
    
    ifPvalues <- file.path(run.name, method, "gpvals.pat.c.%s")
    ifZscores <- file.path(run.name, method, "zmat.%s")	
    ifIndepAnnInputRegion <- file.path(run.name, method, "indep.ann.data.%s")
    ifIndepAbsPos <- file.path(run.name, "data", "abs.start.indep")
    
    ofTable1 <- file.path(run.name, "top.indep.features", 
            paste("TopIndepFeatures%s", method, adjust.method, format(significance, scientific = TRUE, digits=1), ".txt", sep=""))
    
    ##CHECK FILE EXISTENCE
    checkFiles(c(run.name, ifPvalues, ifZscores, ifIndepAnnInputRegion, ifIndepAbsPos))
    
    abs.start.indep.whole <- dget(ifIndepAbsPos)
    
    if(!is.null(input.region.indep)){
        if(length(input.region.indep) != 1)
            stop("Only one input region for the independent data is allowed. The input is now: ", paste(input.region.indep, collapse=", "))
        input.region.indep <- unlist(sapply(input.region.indep, predefinedRegions, USE.NAMES=FALSE))
        region.indep <- getGenomicRegion(input.region.indep)		
        indices.indep <- (abs.start.indep.whole >= region.indep$absolute.start) & (abs.start.indep.whole <= region.indep$absolute.end)
    }
    
    input.regions <- unlist(sapply(input.regions, predefinedRegions, USE.NAMES=FALSE))
    
    rList <- list()
    for(input.region in input.regions)
    {
        region.dep <- getGenomicRegion(input.region)		
        
        cat("... input region:", input.region, "\n")
        
        #find region independent if given
        if(is.null(input.region.indep))
            indices.indep <- (abs.start.indep.whole >= region.dep$absolute.start) & (abs.start.indep.whole <= region.dep$absolute.end)
        
        abs.start.pos <- abs.start.indep.whole[indices.indep]
        
        # Retrieve the p-values based on the provided adjust method. The p-values 
        # are ordered by the position of their corresponding aCGH feature on the
        # chromosome.
        raw.pvals <- dget(sprintf(ifPvalues, input.region))
        
        # Correct the p-values for multiple testing.
        p.values <- p.adjust(raw.pvals, method=adjust.method)
        
        annotation <- dget(sprintf(ifIndepAnnInputRegion, input.region))		
        
        # Retrieve the z-scores.
        load(sprintf(ifZscores, input.region))
       
        z.scores <- z.scores[p.values <= significance , , drop=FALSE]
        
        if(nrow(z.scores) == 0)
            next
        
        mean.influences <- apply(z.scores, 2, mean, na.rm=TRUE)
        
        nan.indices <- is.nan(mean.influences)
        
        mean.influences <- mean.influences[!nan.indices]
        annotation <- annotation[!nan.indices, , drop=FALSE]        
        abs.start.pos <- abs.start.pos[!nan.indices]
        
        
        if(is.null(z.threshold))
            z.threshold <- range(mean.influences)
        
        id.extreme.influences <- mean.influences <=  z.threshold[1] | mean.influences >=  z.threshold[2]
        
        mean.influences <- mean.influences[id.extreme.influences]
        
        abs.start.pos <- abs.start.pos[id.extreme.influences]
        
        annotation <- annotation[id.extreme.influences, , drop=FALSE]
        
        # Sort table by mean influences.
        table <- data.frame(`mean Z-scores`=mean.influences,  `absolute start position`=abs.start.pos, annotation,
                            row.names=NULL)
        
        table <- table[order(table[,1], decreasing=decreasing), ]
        table[,1] <- round(table[,1], digits=2)		
        
        #remove NA's in case of window or overlap
        table <- na.omit(table)        
        
        # Write file to disk.
        write.table(table, file=sprintf(ofTable1, input.region), row.names=FALSE, quote=FALSE,sep="\t")
        
        # Add to return value.
        cat("... tabulated results stored with file name: \n", sprintf(ofTable1, input.region), "\n", sep="")	
        #table is input for make.indep.region()
        rList[[input.region]] <- table
    }
    invisible(rList)
}

