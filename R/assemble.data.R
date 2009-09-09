###############################################################################
# function: checkfiles()
# description: general file existence checker. 
###############################################################################
checkFiles <- function(files)
{
    #maybe options(warnings=FALSE)
    exists <- logical(length(files))
    for(i in 1:length(files))
    {
        if(length(grep("%s", files[i])) > 0){
            pattern <- sub("%s", "", basename(files[i]))
            listFiles <- list.files(path=dirname(files[i]), pattern=pattern)
            exists[i] <- length(listFiles) > 0			
        } else {				
            exists[i] <- file.exists(files[i])
        }				
    }
    
    if(all(!exists))
        stop("Did you run the 'assemble.data' or is the 'run.name' correct?")
    else if(any(!exists) & length(files) > 1) {
        if(grepl("zmat", files[which(!exists)]))
            cat("Probably you will need the Z-scores for this function to run, run the 'integrated.analysis()' with 'zscores=TRUE'!\n")
        stop("File or directory doesn't exists! Filename: ", paste(files[which(!exists)], collapse=", "), "\n", sep="")
    }
    exists
}

###############################################################################
# function: validColumn()
# description: checks if a column either character vector, indices or ligical vector
#              is valid according to given column names.    
###############################################################################
validColumn <- function(x, colNames)
{
    valid <- logical(length(x))
    for(i in 1:length(x))
    {
        xi <- type.convert(as.character(x[[i]]), as.is=TRUE)	
        if(is.character(xi) & all(xi %in% colNames)) {		
            xi <- which(colNames %in% xi)	
            valid[i] <- TRUE
        }	
        if(is.integer(xi) & all(xi > 0 & xi <= length(colNames))) {					
            valid[i] <- TRUE	
        }	
        if(is.logical(xi) & (length(xi) == length(colNames)))	{	
            xi <- which(xi == TRUE)						
            valid[i] <- TRUE
        }	
    }	
    if(any(!valid))
        stop("Invalid column(s): ", paste(x[which(!valid)], collapse=", "))
    xi
}

###############################################################################
# function: assemble.data()
# description: assembles the data necessary to perform the integrate analysis on
###############################################################################
assemble.data <- function(dep.data, 
        indep.data,		
        dep.id="ID", 
        dep.chr="CHROMOSOME", 
        dep.pos="STARTPOS",
        dep.ann=NULL,		
        dep.symb,  
        indep.id="ID", 
        indep.chr="CHROMOSOME", 
        indep.pos="STARTPOS",
        indep.ann=NULL,
        indep.symb,
        overwrite=FALSE, 
        run.name="analysis_results")
{
    cat("Assembling data ...\n")
    
    rescale <- 1e9    
    data("chrom.table", package="SIM")
    
    #save
    ofDepData <- file.path(run.name, "data", "dep.data")
    ofIndepData <- file.path(run.name, "data", "indep.data")
    
    #dput
    ofIndepAbsStart <- file.path(run.name, "data", "abs.start.indep")
    ofDepAbsStart <- file.path(run.name, "data", "abs.start.dep")
    
    ofIndepPos <- file.path(run.name, "data", "indep.pos")
    ofDepPos <- file.path(run.name, "data", "dep.pos")
    
    ofIndepChr <- file.path(run.name, "data", "indep.chr.data")
    ofDepChr <- file.path(run.name, "data", "dep.chr.data")
    
    ofIndepSym <- file.path(run.name, "data", "indep.sym.data")
    ofDepSym <- file.path(run.name, "data", "dep.sym.data")
    
    ofDepAnn <- file.path(run.name, "data", "dep.ann.data")
    ofIndepAnn <- file.path(run.name, "data", "indep.ann.data")
    
    ofDepId <- file.path(run.name, "data", "dep.id.data")
    ofIndepId <- file.path(run.name, "data", "indep.id.data")
    
    #check the input	
    if(!is.data.frame(dep.data) | !is.data.frame(indep.data))		
        stop("The dependent and independent datasets should be of type 'data.frame'.")
    
    if(file.exists(run.name) & !overwrite)		
        stop("Directory ", run.name, " already exist but argument 'overwrite=FALSE'.")
    
    ##check column names	
    validColumn(list(dep.pos, dep.chr, dep.ann), colnames(dep.data))	
    validColumn(list(indep.pos, indep.chr, indep.ann), colnames(indep.data))	
    
    dep.pos <- validColumn(dep.pos, colnames(dep.data))
    dep.pos.data <- dep.data[, dep.pos]
    
    indep.pos <- validColumn(indep.pos, colnames(indep.data))
    indep.pos.data <- indep.data[, indep.pos]
    
    dep.chr <- validColumn(dep.chr, colnames(dep.data))	
    dep.chr.data <- dep.data[, dep.chr]
    
    indep.chr <- validColumn(indep.chr, colnames(indep.data))
    indep.chr.data <- indep.data[, indep.chr]
    
    if(is.null(dep.ann))
        stop("Column name for annotation of the dependent data is missing.")
    dep.ann <- validColumn(list(dep.ann), colnames(dep.data))
    
    if(is.null(indep.ann))
        stop("Column name for annotation independent data is missing")
    indep.ann <- validColumn(list(indep.ann), colnames(indep.data))
    
    ##check column content
    if(!is.numeric(dep.pos.data) | !any(is.finite(dep.pos.data)) & any(is.na(dep.pos.data)))
        stop("The position column of dependent data contains one or more invalid positions.")
    
    if(!is.numeric(indep.pos.data) | !any(is.finite(indep.pos.data)) & any(is.na(indep.pos.data)))
        stop("The position column of independent data contains one or more invalid positions.")
    
    #e.g. no information for mt in chrom.table
    if(!all(levels(factor(toupper(dep.chr.data))) %in% levels(chrom.table$chr)))
        stop("The chromosome column of dependent data contains one or more invalid chromosome identifiers.")
    
    if(!all(levels(factor(toupper(indep.chr.data))) %in% levels(chrom.table$chr)))
        stop("The chromosome column of independent data contains one or more invalid chromosome identifiers.")
    
    #get the right levels for the chromosomes
    dep.chr.data <- factor(dep.chr.data, levels=levels(chrom.table$chr))
    indep.chr.data <- factor(indep.chr.data, levels=levels(chrom.table$chr))
    
    #get the absolute start column for the dependent and independent data by using the
    #user-defined chromosome and basepair position columns	
    abs.start.indep <- as.numeric(indep.pos.data) + as.integer(indep.chr.data) * rescale 	
    abs.start.dep <- as.numeric(dep.pos.data) + as.integer(dep.chr.data) * rescale 	
    
    #order all columns by absulote start and save them in a object with dput() function
    abs.start.indep.ordered <- abs.start.indep[order(abs.start.indep)]
    abs.start.dep.ordered <- abs.start.dep[order(abs.start.dep)]
    
    dep.data <- dep.data[order(abs.start.dep), ]
    indep.data <- indep.data[order(abs.start.indep), ]
    
    dep.pos.data <- dep.pos.data[order(abs.start.dep)]
    indep.pos.data <- indep.pos.data[order(abs.start.indep)]
    
    dep.chr.data <- dep.chr.data[order(abs.start.dep)]
    indep.chr.data <- indep.chr.data[order(abs.start.indep)]
    
    #Why NA's when missing or empty?
    dep.id.data <- NA
    if(!missing(dep.id))
    {
        dep.id <- validColumn(dep.id, colnames(dep.data))
        dep.id.data <- dep.data[, dep.id]
    }
    dep.id.data[dep.id.data == " "] <- NA
    
    indep.id.data <- NA
    if(!missing(indep.id))
    {
        indep.id <- validColumn(indep.id, colnames(indep.data))
        indep.id.data <- indep.data[, indep.id]
    }
    indep.id.data[indep.id.data == " "] <- NA
    
    #generate the necessary folders to put in the data, the results of the analysis and the outputs of the functions
    showWarnings <- ifelse(overwrite, FALSE, TRUE)
    dir.create(run.name, showWarnings=showWarnings)
    dir.create(file.path(run.name, "data"), showWarnings=showWarnings)
    dir.create(file.path(run.name, "pvalue.plots"), showWarnings=showWarnings)
    
    dir.create(file.path(run.name, "heatmap.zscores"), showWarnings=showWarnings)
    dir.create(file.path(run.name, "top.indep.features"), showWarnings=showWarnings)
    dir.create(file.path(run.name, "top.dep.features"), showWarnings=showWarnings)
    
    dput(abs.start.indep.ordered, file=ofIndepAbsStart)
    dput(abs.start.dep.ordered, file=ofDepAbsStart)
    
    dput(indep.pos.data, file=ofIndepPos)
    dput(dep.pos.data, file=ofDepPos)
    
    dput(indep.chr.data, file=ofIndepChr)
    dput(dep.chr.data, file=ofDepChr)
    
    dput(as.vector(dep.id.data), file=ofDepId)
    dput(as.vector(indep.id.data), file=ofIndepId)
    
    dep.ann <- dep.data[, dep.ann]	
    dput(as.matrix(dep.ann), file=ofDepAnn)
    
    indep.ann <- indep.data[, indep.ann]		
    dput(as.matrix(indep.ann), file=ofIndepAnn)
    
    if(!missing(dep.symb))
    {
        dep.symb <- validColumn(dep.symb, colnames(dep.data))
        dep.symb.data <- dep.data[, dep.symb]		
        dput(as.vector(dep.symb.data), file=ofDepSym)
        colnames(dep.data)[dep.symb] <- "SYMBOL"
    }	
    
    if(!missing(indep.symb))
    {
        indep.symb <- validColumn(indep.symb, colnames(indep.data))
        indep.symb.data <- indep.data[, indep.symb]		
        dput(as.vector(indep.symb.data), file=ofIndepSym)
        colnames(indep.data)[indep.symb] <- "SYMBOL"
    }
    
    #set the column names by setting an universal column name for each necessary column
    colnames(dep.data)[c(dep.pos, dep.chr, dep.id)] <- c("STARTPOS", "CHROMOSOME", "ID")	
    colnames(indep.data)[c(indep.pos, indep.chr, indep.id)] <- c("STARTPOS", "CHROMOSOME", "ID")
    
    save(dep.data, file=ofDepData)
    save(indep.data, file=ofIndepData)
    
    cat("... assembled dependent data: dim(", paste(dim(dep.data), collapse=", ", sep=""), ").\n", sep="")
    cat("... assembled independent data: dim(", paste(dim(indep.data), collapse=", ", sep=""), ").\n", sep="")
}

