tabulate.top.dep.features <- function(input.regions="all chrs", 
        adjust.method="BY", 
        method=c("full", "smooth", "window", "overlap"), 
        significance=1,		
        run.name="analysis_results")
{ 
    cat("Tabulate results on dependent data ...\n")
    
    ##CHECK INPUT ARGUMENTS	
    method <- match.arg(method)	
    
    #check significance
    if(!is.numeric(significance) | (significance < 0) | (significance > 1))
        stop("'significance' parameter should be a numeric value between 0 and 1! The input is: ", significance)
    
    ifPvalues <- file.path(run.name, method, "gpvals.pat.c.%s")
    ifZscores <- file.path(run.name, method, "zmat.%s") 
    ifDepAnnInputRegion <- file.path(run.name, method, "dep.ann.data.%s")
    ifIndepAnnInputRegion <- file.path(run.name, method, "indep.ann.data.%s")   
    ifDepAbsPos <- file.path(run.name, "data", "abs.start.dep")
    
    ofTable1 <- file.path(run.name, "top.dep.features", 
            paste("TopDepFeatures%s", method, adjust.method, format(significance, scientific = TRUE, digits=1), ".txt", sep=""))
    
    ##CHECK FILE EXISTENCE
    checkFiles(c(run.name, ifPvalues, ifZscores, ifDepAnnInputRegion, ifDepAbsPos))
    
    #get the absolute start	
    abs.start.dep.whole <- dget(ifDepAbsPos)
    
    input.regions <- unlist(sapply(input.regions, predefinedRegions, USE.NAMES=FALSE))
    
    rList <- list()
    for(input.region in input.regions)
    {
        region.dep <- getGenomicRegion(input.region)		
        
        cat("... input region:", input.region, "\n")
        
        # Retrieve the p-values based on the provided adjust method. The p-values 
        # are ordered by the position of their corresponding dep feature on the
        # chromosome.
        
        #find region independent		
        indices.dep <- (abs.start.dep.whole >= region.dep$absolute.start) & (abs.start.dep.whole <= region.dep$absolute.end)
        abs.start.pos <- abs.start.dep.whole[indices.dep]
        
        raw.pvals <- dget(sprintf(ifPvalues, input.region))				
        p.values <- p.adjust(raw.pvals, method=adjust.method)
        
        annotationDep <- dget(sprintf(ifDepAnnInputRegion, input.region))
        
        annotationIndep <- dget(sprintf(ifIndepAnnInputRegion, input.region))   
        
        load(sprintf(ifZscores, input.region))
        
        #if complete row contains NA remove from data
        na.indices <- apply(z.scores, 1, function(x) sum(is.na(x)) == length(x))
        
        z.scores <- z.scores[!na.indices, , drop=FALSE]
        p.values <- p.values[!na.indices]
        annotationDep <- annotationDep[!na.indices,, drop=FALSE]
        
        abs.start.pos <- abs.start.pos[!na.indices] 
        
        #subset significant P-values
        id.p.values <- p.values <= significance
        
        z.scores <- z.scores[id.p.values , , drop=FALSE]
        
        if(nrow(z.scores) == 0)
            next
        
        annotationDep <- annotationDep[id.p.values, , drop=FALSE]
        abs.start.pos <- abs.start.pos[id.p.values]		
        p.values <- p.values[id.p.values]		
        
        extreme.influences <- apply(z.scores, 1, function(x) x[which.max(abs(x))]) #get all the extreme contributions
        
        extreme.influences.ids <- apply(z.scores, 1, function(x) which.max(abs(x))) #get all the extreme contributions
        
        annotationIndep <- annotationIndep[extreme.influences.ids, ]
        
        # Sort table by the p-values.
        table <- data.frame(`P-values`=p.values, `extreme influences`=extreme.influences,							
                `absolute start position`=abs.start.pos, annotationDep, annotationIndep, row.names=NULL)
        
        table <- table[order(table[,1]), ]
        table[,1] <- signif(table[,1], digits=2)	
        table[,2] <- round(table[,2], digits=2)	
        
        #remove NA's in case of window or overlap
        table <- na.omit(table)        
        
        write.table(table[], file=sprintf(ofTable1, input.region), row.names=FALSE, quote=FALSE, sep="\t")
        
        cat("... tabulated P-values stored with file name: \n", sprintf(ofTable1, input.region), "\n", sep="")	
        #table is input for make.dep.region()
        rList[[input.region]] <- table
    }
    invisible(rList)
}

