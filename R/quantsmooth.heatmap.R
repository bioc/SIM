`quantsmooth.heatmap` <-
function(input.region, dep.data, lambda, run.name)  {
load(sprintf("%s/data/dep.data.only",run.name))

for (input.region.name in input.region) {
    if (substr(input.region.name, 1,3) == "chr")
    {
       chrom_full = strsplit(input.region.name, "_")[[1]][1]
       start_end = strsplit(input.region.name, "_")[[1]][2]
       chrom = as.integer(strsplit(chrom_full, "chr")[[1]][2])
    
       start = as.integer(strsplit(start_end, "-")[[1]][1])
       end = as.integer(strsplit(start_end, "-")[[1]][2])
    
       abs.start = (chrom*1000000000)+start
       abs.end   = (chrom*1000000000)+end
    }
    else if (!substr(input.region.name, 1,3) == "chr")
    {
       options(warn=-1)
       if(is.na(as.integer(input.region.name))) {
          abs.start = chrom.table[input.region.name, "Abs.start"]
          abs.end   = chrom.table[input.region.name, "Abs.end"]
       }

       if(!is.na(as.integer(input.region.name))) {
          chr.pos <- chrom.table[input.region.name == chrom.table$name,]

          if(nrow(chr.pos)==2){
                   abs.start = chr.pos[1, "Abs.start"]
                   abs.end = chr.pos[2, "Abs.end"]
		}
		if(nrow(chr.pos)==1){
                   abs.start = chr.pos[1, "Abs.start"]
                   abs.end = chr.pos[1, "Abs.end"]
		}
       }
       options(warn=1)
    }
    
    abs.start.dep <- dget(sprintf("%s/data/abs.start.dep", 
            run.name))
        indices.dep = (abs.start.dep >= abs.start) & (abs.start.dep <= abs.end)

       if (sum(indices.dep) > 2) {
    
      dep.matrix = as.matrix(dep.data.only[indices.dep, ])
      
     }
	dep.pos.data <- dget(sprintf("%s/data/dep.pos.data", run.name))
	dep.pos <- dep.pos.data[indices.dep]
      dep.pos2 <- seq(0,length(dep.pos)-1,1)
	my.segment = min(100,length(dep.pos2))

	smooth.matrix <- dep.matrix
	for (sample in 1:ncol(dep.matrix))
	{           
	smooth.matrix [,sample] <- quantsmooth(dep.matrix[,sample],smooth.lambda<- lambda, segment = my.segment)
} 

	zmax = max(0, abs(smooth.matrix))      

	xlim <- c(zmax, -zmax)
	cols <- 1:ncol(dep.matrix) + 1
	plot(smooth.matrix[,1],dep.pos2, xlim = xlim, ylim = c(0, max(dep.pos2)), xaxs = 'i', yaxs = 'i', xlab = 'intensity', type = 'l')
	for (sample in 1:ncol(dep.matrix))
	{           
       lines(smooth.matrix [,sample],dep.pos2,col = cols[sample])
	   
     }
     segments(0, 0, 0,max(dep.pos2), col="gray", lty="dotted")
     }   
 }

