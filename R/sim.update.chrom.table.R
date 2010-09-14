sim.update.chrom.table <- function(db="homo_sapiens_core_40_36b") {                                                     
    warning("Be sure RMySQL is loaded!")
    
    ensemblDb <- dbConnect(MySQL(), user="anonymous", db=db, host="ensembldb.ensembl.org")
    
    query <- paste("select name, seq_region_start, seq_region_end, band, stain",  
            "from seq_region s, karyotype k", 
            "where s.seq_region_id=k.seq_region_id;", sep=" ")
    
    table <- dbGetQuery(ensemblDb, query)
    
    #relevel chromosome column first numeric, X, Y
    chrom.table <- data.frame(chr=factor(table$name, level=c(1:22, "X","Y")), 
            arm=substr(table[,"band"], 1, 1), 
            band=gsub("p|q", "", table[,"band"]), 
            start=table$seq_region_start, 
            end=table$seq_region_end, 
            stain=table$stain)
    
    #now order according to chromosome number and start position band
    chrom.table <- chrom.table[order(chrom.table[, "chr"], chrom.table[, "start"]), ]
    
    chrom.table                                                              
}                                                                                




