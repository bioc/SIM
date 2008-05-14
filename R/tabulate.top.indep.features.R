`tabulate.top.indep.features` <-
function(input.regions = "all chrs", adjust.method=c("BY", "BH", "raw"), 
         significance=0.2, sort.order='positive', run.name = NULL)
{ 
  data(chrom.table)
  miss_arms = c("13p", "14p", "15p", "21p", "22p")
  ret = list()
  if(is.null(run.name))
  {
	run.name = "analysis_results"
  }

  if (!file.exists(run.name)) {
    msg = paste("The directory", run.name, "does not exist, please insert a correct run.name.")
    stop(msg)
  }
  if(!file.exists(sprintf("%s/data/indep.data",run.name)))
   {
	msg = "please first run the read.data() function to read in the data necessary for this function"
	stop(msg)
   }

  if(!length(dir(path = sprintf("%s/intermediate.data/", run.name), pattern = "zmat.")) >0)
  {
     stop("no zscores are found. Please run the integrated.analysis() function and set zscores=T")
  }
  load(sprintf("%s/data/indep.data",run.name))

  options(warn = -1)
  if(input.regions == "all chrs"){  input.regions = c(1:24)  }  
  else if (input.regions == "all chrs auto") {input.regions = c(1:22)}  
  else if (input.regions == "all arms") {
        input.regions = rownames(chrom.table)
	input.regions <- input.regions[!(input.regions %in% miss_arms)]
    }
    else if (input.regions == "all arms auto") {
        input.regions = rownames(chrom.table)[1:44]
	input.regions <- input.regions[!(input.regions %in% miss_arms)]
    }
  options(warn = 1)

  for (input.region.name in input.regions) {
   if(input.region.name == "whole genome auto")
   {
      abs.start= min(chrom.table[1:44,]$Abs.start)
      abs.end= max(chrom.table[1:44,]$Abs.end)
   }

   #same for the whole genome
   else if (input.region.name == "whole genome")
   {
      abs.start= min(chrom.table$Abs.start)
      abs.end= max(chrom.table$Abs.end)
   }

   else if (substr(input.region.name, 1,3) == "chr")
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
    

    
    # Retrieve the p-values based on the provided adjust method. The p-values 
    # are ordered by the position of their corresponding aCGH feature on the
    # chromosome.
    options(warn = -1)
    raw.pvals = try( get.pvals(input.region.name,run.name),silent=TRUE )
    options(warn = 1)
    if ( class(raw.pvals) != "try-error" ) {
      adjust.method = match.arg(adjust.method, c("BY", "BH", "raw"))
      
      if (adjust.method == "BY")
        p.values =  adj.byfdr(raw.pvals)
      
      else if (adjust.method == "BH")
        p.values = adj.bhfdr(raw.pvals)
      
      else if (adjust.method == "raw")
        p.values = raw.pvals
        
      
      # Retrieve the z-scores.
      z.scores = get.indep.zscores(input.region.name, run.name)
    
      # Select only those rows in the matrix of z-scores that correspond to
      # significant aCGH features.
      z.scores = z.scores[p.values <= significance, ]

      if ( !is.null(nrow(z.scores)) ) {
        if (nrow(z.scores) > 0) {
          writeln( sprintf("Creating table for %s", input.region.name) )

          # Compute the mean influence of each column (i.e. gene or sample) on
          # the p-values of the significant aCGH features.
          mean.influences = apply(z.scores, 2, mean)
      
          # Select only those rows (aCGH features) that correspond to the 
          # current chromosome and order them according to their absolute start.
          abs.start.indep <- dget(sprintf("%s/data/abs.start.indep",run.name))
   
   	    indices = (abs.start.indep >= abs.start) & (abs.start.indep <= abs.end)

          # Create a data.frame of the p-values.
          selected.columns = data.frame(mean.influences)
          
   	    include.columns = dget(sprintf("%s/data/ann.indep.data",run.name))
     	    include.columns.region = include.columns[indices, ]
          selected.columns = data.frame(selected.columns, include.columns.region)

          # Set the row- and column names.
          colnames(selected.columns) = c("Mean z-score", colnames(include.columns))
                      
          # Sort selected.columns by the p-values.
          if (sort.order == "positive")
            indices = order(selected.columns$"Mean z-score", decreasing=TRUE)
          else if (sort.order == "negative")
            indices = order(selected.columns$"Mean z-score")
	  else{
	  msg = paste("The inserted sort.order must be negative or positive, you've inserted", sort.order)
          stop(msg)
          }
            
          selected.columns.sorted = as.data.frame( selected.columns[indices, ] )
          colnames(selected.columns.sorted) = colnames(selected.columns)
                
          # Write file to disk.
          filename = sprintf("%s/top.indep.features/top.indep.features.%s.txt",run.name, input.region.name)
          write.table(selected.columns.sorted, file=filename, row.names=FALSE, quote = FALSE,sep = "\t")
          
          # Add to return value.
          ret[[input.region.name]] = selected.columns.sorted
          writeln()
          writeln("**************************************************")
          writeln("The results are stored in the following file:")
          writeln(filename)
          writeln("**************************************************")
          writeln()
          
        }
        
    }
  }
  else {
      writeln("Will not generate table for '%s': it has no features.", input.region.name)
      writeln()
    }
  invisible(ret)    
}
}

