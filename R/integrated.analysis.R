###############################################################################
# function: weightMatrix()
# description: constructs and weight-matrix for the withinWindow-function 
###############################################################################
weightMatrix <- function(x, n)
{
    x <- as.matrix(x)	
    if(nrow(x) == 1 & ncol(x) == 1){
        x <- matrix(rep(x, 2*n), ncol=2)
    }else if(ncol(x) == 1 & nrow(x) == 2){
        x <- matrix(rep(x, each=n), ncol=2)
    }else if(ncol(x) == 1 & nrow(x) == n){
        x <- matrix(rep(x, 2), ncol=2)
    }
    x
}

###############################################################################
# function: withinWindow()
# description: given two vectors containing base pair position x, y and a weight-matrix, w
#              the positions in x that are in [y[i]-w[i], y[i]-w[i]] is returned. 
###############################################################################
withinWindow <- function(x, y, w, sort=TRUE)
{
    ##different window options are:
    #overlap:                                                           w=0	 
    #fixed window size for each y but left and right are equal:         w > 0  
    #fixed window size for each y but left and right are different:     w=c(10, 100)
    #different window size for each y but left and right are equal:     nrow(w)=length(y)
    #different window size for each y but left and right are different: dim(w)=c(length(y), 2)
    
    w <- weightMatrix(w, n=length(y))
    
    #maybe other checks
    if(sort)
    {
        x <- sort(x)
        y <- sort(y)	
    }
    
    w <- .C("withinWindow", x=as.vector(x, mode="double"), 
            y=as.vector(y, mode="double"),
            nx=as.integer(length(x)),			 
            ny=as.integer(length(y)),
            w=as.vector(w, mode="double"), PACKAGE="SIM")$w
    w <- matrix(w, nrow=length(y), ncol=2)
    rownames(w) <- y
    w
}

###############################################################################
# function: predefiendRegions()
# description: transforms predefined input region to chromosomal regions
#              some chromosomal region are less interesting and are excluded
###############################################################################
predefinedRegions <- function(x, excluded.arms=c("13p", "14p", "15p", "21p", "22p"))
{ 
    x <- as.character(x)
    #convert UCSC format to default
    x <- gsub("^chr", "", x)
    x <- gsub(":", " ", x)	
    switch(x,
            `whole genome`=c("1p", "Yq"),
            `whole genome auto`=c("1p", "22q"),
            `all chrs`=c(1:22, "X", "Y"),
            `all chrs auto`=c(1:22),
            `all arms`=setdiff(c(paste(c(1:22, "X", "Y"), "p", sep=""), paste(c(1:22, "X", "Y"), "q", sep="")) , excluded.arms),
            `all arms auto`=setdiff(c(paste(c(1:22), "p", sep=""), paste(c(1:22), "q", sep="")), excluded.arms),
            x)
}

###############################################################################
# function: convertGenomicRegion()
# description: converts the a genomic region to start and end base pair position
#              using a chromosome table.
###############################################################################
convertGenomicRegion <- function(...) 
{
    data("chrom.table", package="SIM")
    n <- nargs(...)
    args <- list(...)
    arguments <- list(chr=unique(chrom.table$chr), arm=unique(chrom.table$arm), band=unique(chrom.table$band))
    Columns <- c("start", "end")
    Rows <- TRUE
    for(i in 1:n)
    {
        args[i] <- match.arg(as.character(toupper(args[i])), toupper(arguments[[i]]))		
        Columns <- c(Columns, names(arguments[i]))
        Rows <- Rows & toupper(chrom.table[,names(arguments[i])]) == toupper(args[i])      
    }
    
    select <- chrom.table[Rows, Columns]	
    
    if(nrow(select) < 1)
        stop("Invalid combination of either ", paste(names(arguments), sep=", "))
    
    select[1, "end"] <- select[nrow(select),"end"]
    select[1,]
}

###############################################################################
# function: Bases()
# description: converts and checks if a given start and end position is valid
###############################################################################
Bases <- function(chrx, top, bottom) 
{
    top <- as.integer(top)	
    bottom <- as.integer(bottom)	
    selected <- convertGenomicRegion(chrx)
    if(nrow(selected) < 1 | selected$start > top | selected$end < bottom)
        stop("Invalid combination of ", paste(chrx, top, bottom, sep=", "))
    selected$start <- top 
    selected$end <- bottom
    selected
}

###############################################################################
# function: getGenomicRegion()
# description: converts userdefined region to start and end position
###############################################################################
getGenomicRegion <- function(input.region, rescale=1e9) 
{	
    input.region <- as.character(input.region)
    
    extractChr <- function(x) {chr <- unlist(strsplit(x, "[^0-9|xy|XY]{1,2}"))[1]; ifelse(nchar(chr) > 0, chr, NA)}
    extractArm <- function(x) {loc <- regexpr("[p|q]", x); ifelse(loc > 0, substr(x, loc, loc), NA)}
    extractBand <- function(x) {loc <- regexpr("[p|q]", x); ifelse(loc > 0 & loc != nchar(x), substr(x, loc+1, nchar(x)), NA)}
    extractBases <- function(x) {loc <- regexpr(" [0-9]+-[0-9]+", x); unlist(ifelse(loc > 0, strsplit(substr(x, loc+1, nchar(x)),"-"), NA))}
    
    chr <- extractChr(input.region)
    arm <- extractArm(input.region) 
    band <- extractBand(input.region)
    bases <- extractBases(input.region)
    region <- NA
    if(!is.na(chr) & !is.na(arm) & !is.na(band))
        region <- convertGenomicRegion(chr, arm, band) 
    else if(!is.na(chr) & !is.na(arm) & is.na(band))
        region <- convertGenomicRegion(chr, arm)
    else if(!is.na(chr) && !is.na(bases))
        region <- Bases(chr, bases[1], bases[2])
    else if(!is.na(chr) & is.na(arm) & is.na(band))
        region <- convertGenomicRegion(chr)	
    else
        stop("Unknown input region: ", input.region)		
    
    #add absolute genomic scale
    region$absolute.start <- as.numeric(region$start) + as.numeric(region$chr)*rescale 
    region$absolute.end <- as.numeric(region$end) + as.numeric(region$chr)*rescale
    region
}

###############################################################################
# function: integrated.analysis()
# description: performes the integrated analysis on two genomic datasets by
#              calling runIA(). 
###############################################################################
integrated.analysis <- function (samples, 
        input.regions="all chrs",
        input.region.indep=NULL,	 
        zscores=FALSE,			
        method=c("full", "smooth", "window", "overlap"),		
        dep.end=1e5, 
        window=c(1e6, 1e6), 
        smooth.lambda=2, 
        adjust=~1,
        run.name="analysis_results", ...)
{
    #convert input regions
    input.regions <- unlist(sapply(input.regions, predefinedRegions, USE.NAMES=FALSE))
    
    cat("Performing integrated analysis on input region(s): ", paste(input.regions, collapse=", ", sep=""), ".\n", sep="")
    
    ##CHECK INPUT ARGUMENTS	
    method <- match.arg(method)
    
    #check input parameters
    if(method == "window" & (!is.numeric(window) | length(window) != 2 | any(window < 0))) 
        stop("'window' parameter is not valid!", window)
    if(method == "smooth" & (!is.numeric(smooth.lambda) | smooth.lambda <= 0)) 
        stop("'smooth.lambda' parameter is not valid!", smooth.lambda)
    
    ##CREATE DIRECTORY FOR METHOD SPECIFIC DATA
    if(!file.exists(file.path(run.name, method)))
        dir.create(file.path(run.name, method))
    
    ##FILES/DIRECTORIES THAT ALREADY EXISTS 
    #load
    ifDepData <- file.path(run.name, "data", "dep.data")
    ifIndepData <- file.path(run.name, "data", "indep.data")
    
    #dget
    ifDepAbsPos <- file.path(run.name, "data", "abs.start.dep")
    ifIndepAbsPos <- file.path(run.name, "data", "abs.start.indep")
    
    ifDepPos <- file.path(run.name, "data", "dep.pos")
    ifIndepPos <- file.path(run.name, "data", "indep.pos")
    
    ifDepAnn <- file.path(run.name, "data", "dep.ann.data")
    ifIndepAnn <- file.path(run.name, "data", "indep.ann.data")
    
    ##CHECK FILE EXISTENCE
    checkFiles(c(ifDepData, ifIndepData, ifDepAbsPos, ifIndepAbsPos, ifDepPos, ifIndepPos))
    
    #save
    ofDepData <- file.path(run.name, "data", "dep.data.only")
    ofIndepData <- file.path(run.name, "data", "indep.data.only")
    
    ##FILES THAT WILL BE CREATED
    #dput %s will be replaced by inputRegion
    ofDepPosInputRegion <- file.path(run.name, method, "start.input.region.%s") 
    ofDepDataInputRegion <- file.path(run.name, method, "data.dep.pat.c.%s")   
    ofIndepDataInputRegion <- file.path(run.name, method, "data.indep.pat.c.%s")	
    ofDepAnnInputRegion <- file.path(run.name, method, "dep.ann.data.%s")
    ofIndepAnnInputRegion <- file.path(run.name, method, "indep.ann.data.%s")
    
    ofZscoresInputRegion <- file.path(run.name, method, "zmat.%s") 
    ofPvaluesInputRegion <- file.path(run.name, method, "gpvals.pat.c.%s")
    
    #get the datasets
    load(ifDepData)
    load(ifIndepData)
    
    ##CHECK INPUT DATA	
    dep.samples <- validColumn(list(samples), colnames(dep.data))
    indep.samples <- validColumn(list(samples), colnames(indep.data))
    
    if(length(dep.samples) != length(indep.samples))
        stop("Selected columns dependent and independent data are not identical!")
    
    #subset the data on selected samples
    oldColnames <- colnames(dep.data)
    indep.data <- indep.data[, indep.samples] 
    dep.data <- dep.data[, dep.samples]
    
    if(!(apply(dep.data, 2, is.numeric) || apply(indep.data, 2, is.numeric)))
        stop("Dependent or independent data contains non-numeric values.")
    
    if(sum(is.na(dep.data)) > 0)	
        stop("NA's are not allowed in the dependent data.\n",
                "Please run the function 'impute.nas.by.surrounding()'", 
                "to replace the NA's by the median of the surrounding features")
    
    save(dep.data, file=ofDepData)
    save(indep.data, file=ofIndepData)
    
    #get some data
    abs.start.indep.whole <- dget(ifIndepAbsPos)
    abs.start.dep.whole <- dget(ifDepAbsPos)
    
    if(!is.null(input.region.indep)){
        if(length(input.region.indep) != 1)
            stop("Only one input region for the independent data is allowed. The input is now: ", paste(input.region.indep, collapse=", "))
        
        input.region.indep <- unlist(sapply(input.region.indep, predefinedRegions, USE.NAMES=FALSE))
        region.indep <- getGenomicRegion(input.region.indep)		
        indices.indep <- (abs.start.indep.whole >= region.indep$absolute.start) & (abs.start.indep.whole <= region.indep$absolute.end)
    }
    
    dep.pos.data <- dget(ifDepPos)
    indep.pos.data <- dget(ifIndepPos)
    
    ann.dep <- dget(ifDepAnn)		
    ann.indep <- dget(ifIndepAnn)	
    
    for(input.region in input.regions)
    {
        region.dep <- getGenomicRegion(input.region)
        
        cat("... input region:", input.region, "\n")
        
        #find region independent		
        indices.dep <- (abs.start.dep.whole >= region.dep$absolute.start) & (abs.start.dep.whole <= region.dep$absolute.end)
        
        #find region independent if given
        if(is.null(input.region.indep))
            indices.indep <- (abs.start.indep.whole >= region.dep$absolute.start) & (abs.start.indep.whole <= region.dep$absolute.end)
        
        if(!(sum(indices.indep) > 0 & sum(indices.dep) > 0)){ 			
            cat("There are no features to run for!\n")
            next
        }
        
        #subset data on selected features
        indep.data.region <- as.matrix(indep.data[indices.indep, ])		
        dep.data.region <- as.matrix(dep.data[indices.dep, ])
        
        x <- abs.start.indep.whole[indices.indep]
        y <- abs.start.dep.whole[indices.dep]
        
        rescale <- 1e9 #this is not so nice
        
        if(method == "full"){
            window <- max(y) #max(x, y)?			
        }else if(method == "overlap"){
            window <- 0			
        }else if(method == "smooth"){
            segment <- min(nrow(dep.data.region), 100)
            dep.data.region <- apply(dep.data.region, 2, quantsmooth, smooth.lambda=smooth.lambda, segment=segment)
            y <- c(1, y[-1])
            window <- matrix(c(rep(0, length(y)), diff(c(y, 2*y[length(y)]))-1), ncol=2)		
        }else if(method == "window"){
            if(!is.numeric(dep.end)){
                dep.end <- validColumn(dep.end, oldColnames)
                endpos <- as.numeric(dget(ifDepAnn)[, dep.end])
                y <- (dep.pos.data[indices.dep] + endpos[indices.dep])/2
                y <- y + as.numeric(region.dep$chr)*rescale	
            }
            else
                y <- y + dep.end		  
        }
        
        subset <- withinWindow(x, y, w=window, sort=FALSE)
        
        #remove empty subsets?
        save(subset, file=file.path(run.name, method, "subset.RData"))
        
        result <- runIA(dep.data.region, indep.data.region, zscores, subset, adjust, ...)	
        
        ##STORE MODIFIED DATA IN METHOD SPECIFIC DIRECTORY
        dput(ann.dep[indices.dep, ], file= sprintf(ofDepAnnInputRegion, input.region))
        dput(ann.indep[indices.indep, ], file= sprintf(ofIndepAnnInputRegion, input.region))
        
        dput(dep.pos.data[indices.dep], file=sprintf(ofDepPosInputRegion, input.region))
        
        dput(dep.data.region, file=sprintf(ofDepDataInputRegion, input.region))
        dput(indep.data.region, file=sprintf(ofIndepDataInputRegion, input.region))
        
        ##STORE ZSCORES AND PVALUES
        dput(result$pvalues, file=sprintf(ofPvaluesInputRegion, input.region))
        if(zscores)
        {
            #dput/dget is very slow
            #dput(result$zscores, file=sprintf(ofZscoresInputRegion, input.region))
            z.scores <- result$zscores
            save(z.scores, file=sprintf(ofZscoresInputRegion, input.region))
        }
    }
    
    #make a log-file for each run
    cat(paste("log file            :", date(), 
                    "\nsamples used      :", paste(samples, collapse=", "),
                    "\ninput region used :", paste(input.regions, collapse=", "),	
                    if(!is.null(input.region.indep)) "\nindependent input region used :", paste(input.region.indep, collapse=", "),	
                    "\nzscores           :", zscores,			  
                    "\nmethod            :", method,
                    "\nadjust            :", paste(adjust, collapse="")),						
            if(method == "window") 
                paste("\ndependent end     :", dep.end,
                        "\nwindow            :", paste(window, collapse=", ")),						 
            file=file.path(run.name, paste("log_file_", format(Sys.time(), "%a_%b_%d_%Hh_%Mm_%Ss_%Y"), ".txt", sep="")))
    
}
