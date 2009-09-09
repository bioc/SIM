`RESOURCERER.annotation.to.ID` <-
        function(data, poslist, col.ID.link=1, col.poslist.link=1){
    #first get only the identifier and the Phy.map column
    poslist2 <- data.frame(poslist[,col.poslist.link],poslist[,"Phy.Map"], poslist[,6])
    
    colnames(poslist2) <- c(colnames(poslist)[col.poslist.link], "Phy.Map", "Symbol")
    
    #merge the identifier column with the data
    dat.merged2 <- merge(poslist2, data, by.x=col.poslist.link, by.y=col.ID.link)
    
    if(nrow(dat.merged2) == 0)
        stop("the inserted columns col.ID.link and col.poslist.link don't have any comparitive probes")
    dat.merged <- dat.merged2[dat.merged2[,2] != "",]
    
    #get only the Phy.Map column
    chrom <- dat.merged[,2]
    chrom <- as.character(chrom)
    #script to get the chromosoom position out of a string:
    #chrom is for example:
    #"chr6 (30,960,320-30,975,908)"
    #so split on the "(" character to get the chromosoom and start+end
    #when splitted the strings are: 
    #"chr6 "     "30,960,320-30,975,908)"
    #first get the chromosoom (="chr6":
    
    anno <- array(data=NA, dim=c(length(chrom),2))
    options(warn=-1)
    for(i in 1:length(chrom)){
        chrom_full=strsplit(chrom[i], "\\(")[[1]]
        #now split between chr en 6:
        chrom_full <- gsub("X", 23, chrom_full)
        chrom_full <- gsub("Y", 24, chrom_full)
        chr=strsplit(chrom_full, "chr")[[1]][2]
        anno[i,1] <- as.integer(strsplit(chr, "_")[[1]][1])
        
        #get the start and end:
        start_end2=strsplit(chrom[i], "\\(")[[1]][2]
        #remove the ")" at the end
        start_end=strsplit(start_end2, "\\)")[[1]]
        
        #now split on the "-" character
        start=strsplit(start_end, "-")[[1]][1]
        
        #remove the comma's
        anno[i,2] <- as.integer(gsub("\\,", "", start))
    }
    options(warn=1)
    
    colnames(anno) <- c("CHROMOSOME", "START_POS")
    
    Symbols <- as.character(dat.merged[,3])
    
    Symbol <- Symbols
    
    for(i in 1:length(Symbols))
        Symbol[i]=strsplit(Symbols[i], "\\;")[[1]][1]
    
    Symbol[is.na(Symbol)] <- ""
    data <- dat.merged[,-c(2:3)]
    
    if(col.ID.link == 1)
        dataset <- data.frame(data[,1], anno, Symbol, data[,3:ncol(data)])
    
    else
        dataset <- data.frame(data[,1:col.ID.link], anno, Symbol, data[,(col.ID.link+2):ncol(data)])
    
    colnames(dataset)[1] <- colnames(dat.merged)[1]
    
    dataset <- dataset[!is.na(dataset$CHROMOSOME),]
    return(dataset)
}

