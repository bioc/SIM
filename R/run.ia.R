###############################################################################
# function: showProgress()
# description: nicely prints seconds into hours, minutes and seconds.
###############################################################################
showProgress <- function(i, n, s, nskip){	
    cat("    ... features", i, "of", n)
    s <- as.numeric(s)*(n-i)/nskip
    h <- s%/%3600
    m <- (s - h*3600)%/%60
    s <- (s - h*3600 - m*60)
    if(h > 0)
        cat(" (", h, "h ", m, "min ", round(s, 0), "s)\n", sep="")
    else if(h == 0 & m > 0)
        cat(" (", m, "min ", round(s, 0), "s)\n", sep="")
    else
        cat(" (", round(s, 0), "s)\n", sep="")
    flush.console() #for windows
}

###############################################################################
# function: runIA()
# description: the interface to globaltest, performs on each dependent and independent
#              feature a globaltest and extract the necessary information from a globaltest-object
#              for further analysis. 
###############################################################################
runIA <- function(Y, X, zscores, subset, adjust, ...) #parameters passed to globaltest e.g. adjust, permutation
{
    associatedZscores <- NULL
    if(zscores)
        associatedZscores <- matrix(NA, nrow=nrow(Y), ncol=nrow(X)) #by default all unknown
    
    pValues <- rep(1, nrow(Y)) #by default all unsignificant
    
    #define number of skips
    nskip <- ifelse(as.logical(floor(nrow(Y)/20)), floor(nrow(Y)/20), 1)
    
    start <- Sys.time()
    for(idx in 1:nrow(Y))
    {
        
        y <- Y[idx, ]		
        
        if(subset[idx, 1] == 0 & subset[idx, 2] == 0) #skip empty subset			
            next
        
        sbst <- do.call(":", as.list(subset[idx,]))
        
        object <- gt(response=y, alternative=X, null=adjust, subsets=sbst, ...)		
        
        #extract information from globaltest object
        
        pValues[idx] <- p.value(object)		
        
        # Determine the z-scores for all genes.
        if (zscores) {			
            
            # get the test function
            test <- function(set) object@functions$test(set, calculateP=FALSE)
            
            # Test covariates  
            leaves <- t(sapply(1:size(object), function(i) test(i)))
            
            # calculate zscores  
            zsc <- (leaves[,"S"]  - leaves[,"ES"]) / leaves[,"sdS"]
            
            # Set z-scores with an "NA" value and all z-scores < 0 to zero.
            zsc[is.na(zsc) | zsc < 0] <- 0
            
            #association
            positive <- object@functions$positive()[sbst]
            
            associatedZscores[idx, sbst] <- zsc * (2*positive - 1) # f(0)=-1; f(1)=1 => f(x)=2*x - 1
        }		
        
        if(idx %% nskip == 0){
            showProgress(idx, nrow(Y), difftime(Sys.time(), start, units="secs"), nskip)
            start <- Sys.time()
        }
        
    }#end for-loop	
    
    if(exists("object"))
    {	
        #get globaltest object information taken from globaltest summary
        df <- object@functions$df()
        nperms <- object@functions$nperms()
        
        cat("    ... summary results of 'integrated.analysis()' on last feature:\n")	
        cat("    ... \"gt.object\" object from package globaltest\n")
        cat("    ... Call:\n")
        cat("    ... ", deparse(object@call), "\n")	
        cat("    ... Model:", object@model, "regression.\n")
        cat("    ... Degrees of freedom:", df[1], "total;", df[2], "null;", df[2], "+", df[3], "alternative.\n")
        cat("    ... Null distibution: ")	
        if (nperms[1]) {
            cat(if(!nperms[2]) "all", nperms[1], if(nperms[2]) "random", "permutations.\n")
        } else {
            cat("asymptotic.\n")
        }
        cat("\n")
    }
    else
        cat("    ... there where no independent features for running the 'globaltest()'!\n")
    list(zscores=associatedZscores, pvalues=pValues)
}

