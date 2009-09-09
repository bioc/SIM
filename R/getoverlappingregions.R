#function that calculates overlap between two inseted tables either by 'union' or 'overlapping'
getoverlappingregions <- function(independent_regions, 
		                            dependent_regions, 
		                            method = c("union", "overlapping"))
{
	overlap.region <- c()
	genes.symb <- c()
	
	if(nrow(independent_regions) == 0 | nrow(dependent_regions) == 0)
		stop("either the independent regions or the dependent regions contains no regions, no overlapping regions can be found")
		
	#load biomaRt to find the genes that are within a found region
	require(biomaRt)
	
	ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
	
	regions.indep <- as.character(independent_regions[,2]) #absolute start positions?
	regions.dep <- as.character(dependent_regions[,2])
	
	
	#UNION
	if(method == "union"){
	
		#go through each found dependent region and cut the start and end
		for(reg in 1:length(regions.dep)){
		
			chrom_full <- strsplit(as.character(regions.dep[reg]), "_")[[1]][1] #get the chromosome			
			start_end <- strsplit(as.character(regions.dep[reg]), "_")[[1]][2]
			start <- as.integer(strsplit(start_end, "-")[[1]][1])
			end <- as.integer(strsplit(start_end, "-")[[1]][2])
			
			#go through each found independent region and cut the start and end
			for(reg2 in 1:length(regions.indep)){
			
				chrom_full2 <- strsplit(as.character(regions.indep[reg2]), "_")[[1]][1]
				start_end2 <- strsplit(as.character(regions.indep[reg2]), "_")[[1]][2]
				start2 <- as.integer(strsplit(start_end2, "-")[[1]][1])
				end2 <- as.integer(strsplit(start_end2, "-")[[1]][2])
				
				#calculate if there is an overlap 
				#if there is, union the overlapping region. 
				#Say independent region = 1-10, dependent region = 5-12. The output is 1-12. 
				
				if(start2 >= start & start2 <= end | end2 >= start & end2 <= end & chrom_full==chrom_full2){
					minreg <- min(c(start, start2, end, end2))
					maxreg <- max(c(start, start2, end, end2))
					full_reg <- sprintf("%s_%s-%s", chrom_full, minreg, maxreg)
					overlap.region <- c(overlap.region, full_reg)
					
					chrom <- as.integer(strsplit(chrom_full, "chr")[[1]][2])
					
					#of the found overlapping region, find the genes that are within that region
					
					out <- getBM(attributes=c("hgnc_symbol"),
							       filters=c("chromosome_name", "start", "end"),
							       values=list(chrom, minreg, maxreg), mart=ensembl)
					
					cvec <- c()
					if(!is.null(out)){
						for(i in 1:nrow(out))
							cvec <- paste(cvec, out[i,], sep="; ")
						genes.symb <- c(genes.symb, cvec)
					}
					else
						genes.symb <- c(genes.symb, NA)
				}
			}
		}
	}
	
	#overlap
	if(method == "overlapping"){
		overlap.region <- c()
		
		#go through each found dependent region and cut the start and end
		
		for(reg in 1:length(regions.dep)){
		
			chrom_full <- strsplit(as.character(regions.dep[reg]), "_")[[1]][1]			
			start_end <- strsplit(as.character(regions.dep[reg]), "_")[[1]][2]
			start <- as.integer(strsplit(start_end, "-")[[1]][1])
			end <- as.integer(strsplit(start_end, "-")[[1]][2])
			
			#go through each found independent region and cut the start and end 
			
			for(reg2 in 1:length(regions.indep)){
			
				chrom_full2 <- strsplit(as.character(regions.indep[reg2]), "_")[[1]][1]
				start_end2 <- strsplit(as.character(regions.indep[reg2]), "_")[[1]][2]
				start2 <- as.integer(strsplit(start_end2, "-")[[1]][1])
				end2 <- as.integer(strsplit(start_end2, "-")[[1]][2])
				
				#calculate if there is an overlap
				#if there is, calculate which part of the regions overlap 
				
				if((pmin(end, end2) - pmax(start, start2)) > 0 & chrom_full == chrom_full2){
					full_reg <- sprintf("%s_%s-%s", chrom_full, pmax(start, start2), pmin(end, end2))
					overlap.region <- c(overlap.region, full_reg)
					
					chrom <- as.integer(strsplit(chrom_full, "chr")[[1]][2])
					
					#of the found overlapping region, find the genes that are within that region
					
					out <- getBM(attributes=c("hgnc_symbol"),
							       filters=c("chromosome_name","start","end"),
							       values=list(chrom,pmax(start, start2),pmin(end, end2)), mart=ensembl)
					
					cvec <- c()
					if(!is.null(out)){
						for(i in 1:nrow(out))
							cvec <- paste(cvec, out[i,], sep="; ")
						genes.symb <- c(genes.symb, cvec)
					}
					else
						genes.symb <- c(genes.symb, NA)
				}
			}
		}
	}
	
	
	overlap.region2 <- data.frame(overlap.region, genes.symb)
	
	if(is.null(overlap.region))
		print("no overlapping regions found")
	else
		return(overlap.region2)
		
}
