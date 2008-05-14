`sim.plot.zscore.heatmap` <-
function(input.regions = "all chrs", significance=0.2,
         z.threshold=3, show.names.indep=FALSE, show.names.dep=FALSE,adjust.method = c("BY", "BH", "raw"), scale = "auto", plot.method = c("none"), Normal.data=if(plot.method =="clac") FALSE, windowsize =5,lambda = 2, subtype = FALSE, acgh.heatmap.scale = "auto", pdf = TRUE, run.name = NULL,...)  {
while(names(dev.cur()) == "pdf") dev.off()

  data(chrom.table)
  zmax = 0

  if(is.null(run.name))
  {
	run.name = "analysis_results"
  }
  if (!file.exists(run.name)) {
    msg = paste("The directory", run.name, "does not exist, please insert a correct run.name.")
    stop(msg)
  }
   if(!file.exists(sprintf("%s/data/dep.data", run.name)))
   {
	msg = "please first run the read.data() function to read in the data necessary for this function"
	stop(msg)
   }

  if(!length(dir(path = sprintf("%s/intermediate.data/", run.name), pattern = "zmat.")) >0)
  {
     stop("no zscores are found. Please run the integrated.analysis() function and set zscores=T")
  }

  if(!is.numeric(windowsize))
	stop("the windowsize should be numeric")

  if(length(scale) == 1)
if(scale != "auto"){ stop("the inserted scale is not correct") }
	
  if(length(scale) == 2 & (scale[1] > scale[2] | scale[1] > 0)) stop("the inserted scale is not correct, the first variable in scale should be lower than the second variable. The first variable should also be higher than 0. ") 

if(length(scale) > 2) stop("the inserted scale contains more than 2 variables")

  if(length(acgh.heatmap.scale) == 1)
if(acgh.heatmap.scale != "auto"){ stop("the inserted acgh.heatmap.scale is not correct") }
	
  if(length(acgh.heatmap.scale) == 2 & (acgh.heatmap.scale[1] > acgh.heatmap.scale[2] | acgh.heatmap.scale[1] > 0)) stop("the inserted acgh.heatmap.scale is not correct, the first variable in acgh.heatmap.scale should be lower than the second variable. The first variable should also be higher than 0.") 

if(length(acgh.heatmap.scale) > 2) stop("the inserted acgh.heatmap.scale contains more than 2 variables")


  load(sprintf("%s/data/dep.data",run.name))
  load(sprintf("%s/data/indep.data",run.name))

  #if asked for a predefined input region, get the input.regions:
  options(warn = -1)
  miss_arms = c("13p", "14p", "15p", "21p", "22p")
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


    
    z.scores = try(get.indep.zscores(input.region.name,run.name), silent=T)

    if ( is.null(attr(z.scores, "class")) ) {
      zmax = max(zmax, abs(z.scores))

    } else {
      writeln("Will not generate heatmap for '%s': it has no features.",
              input.region.name)
      writeln()

    }
 }
  for (input.region.name in input.regions) {

    writeln("Running for input.region: %s", input.region.name)
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
    abs.start.indep <- dget(sprintf("%s/data/abs.start.indep",run.name))
    abs.start.dep <- dget(sprintf("%s/data/abs.start.dep",run.name))

    indices.indep = (abs.start.indep >= abs.start) & (abs.start.indep <= abs.end)
    indices.dep = (abs.start.dep >= abs.start) & (abs.start.dep <= abs.end)
    if (sum(indices.dep) > 2 & sum(indices.indep) > 2) {
      dep.data.input.region= dep.data[indices.dep, ]
      indep.data.input.region = indep.data[indices.indep, ]

      # Load the zscores for the genes from disk. Each row corresponds to an
      # aCGH feature. The columns correspond to the genes.
      z.scores = get.indep.zscores(input.region.name,run.name)
     
     if (nrow(dep.data.input.region) != nrow(z.scores))
        stop("Dimensions of dep.data.input.region and the z-score matrix do not match!")

      # Load the p-values that indicate the significance of the association
      # between *all* genes and the aCGH feature.
      raw.pvals = get.pvals(input.region.name,run.name)
      # Correct the p-values for multiple testing.
      adjust.method = match.arg(adjust.method, c("BY", "BH", "raw"))

      if (adjust.method == "BY")
        p.values =  adj.byfdr(raw.pvals)

      else if (adjust.method == "BH")
        p.values = adj.bhfdr(raw.pvals)

      else if (adjust.method == "raw")
        p.values = raw.pvals

      #
      # FALSEirst, determine which genes *on average* have the largest impact on the
      # p-values. Only features with a significant p-value are used to compute the
      # the mean z-scores.
      #

      z.scores.sig = z.scores[p.values <= significance, ]

      # Compute the mean influence each gene has on the p-value of an aCGH feature.

      if ( is.null(dim(z.scores.sig)) )
        z.scores.mean = z.scores.sig

      else
        z.scores.mean = apply(z.scores.sig, 2, mean)

      #
      # Next, create the heatmap.
      #

      # Create a new vector that has as many entries as there are genes. Set the
      # entries that are in the previously computed list of 10 lowest scoring
      # genes to "1" and those that are in the list of 10 highest to "2".

      zscores.vec = rep(1, length(z.scores.mean))
      zscores.vec[ z.scores.mean < -z.threshold ] = 0
      zscores.vec[ z.scores.mean >  z.threshold ] = 2

      # Create a list of colors from the vector created above. This causes the
      # genes in the top 10 lists to be differentially colored.

      levels = 0:2
      colors = c("#ff0000", "#8f8f8f", "#00ff00")
      vec.col.top = create.color.vector(zscores.vec, levels, colors)[[1]]

	if(significance > 0.20){
	
p.values.discrete <- p.values
p.values.discrete [p.values.discrete  > significance] <- 0
p.values.discrete [(p.values >0.2) & (p.values <= significance)] <- 1 
p.values.discrete[(p.values > 0.05) & (p.values <= 0.20)] <- 2
p.values.discrete[p.values <= 0.05] <- 3

      levels = 0:3
      colors = c( "white", "#9f9f9f", "#0088ff", "#00ffff")
      color.list    = create.color.vector(p.values.discrete, levels, colors)
pval.legend = c("p > 0.20", paste("0.20 > p <= ", significance), "0.05 < p <= 0.20", "p <= 0.05") 

}
else{
      p.values.discrete = 0 * (p.values > 0.20) +
                              1 * ((p.values > 0.05) & (p.values <= 0.20)) +
                              2 * (p.values <= 0.05)

      levels = 0:2
      colors = c("white", "#0088ff", "#00ffff")
      color.list    = create.color.vector(p.values.discrete, levels, colors)
pval.legend = c("p > 0.20", "0.05 < p <= 0.20", "p <= 0.05")
}

      # Create a color-vector of the discrete p-values yielded by the globaltest
      # (adjusted for multiple testing).

      vec.col.pvals = color.list[[1]]
      vec.col       = color.list[[2]]

      # Clear/set the column names depending on the number of samples..
      rownames(z.scores) = rep("", nrow(z.scores))
      colnames(z.scores) = rep("", ncol(z.scores))

      if (show.names.dep) {
	  if(!file.exists(sprintf("%s/data/dep.id.data",run.name)))
   	  {
		 msg = paste("please first run the read.data() function to read in the ID's of the data necessary for this function")
    	 	 stop(msg)
   	  }

	  dep.id_all <- dget(sprintf("%s/data/dep.id.data",run.name))
	  dep.id <- dep.id_all[indices.dep]
        new_names = rep("", nrow(z.scores))
        if(!file.exists(sprintf("%s/data/dep.symb.data",run.name)))
	  {
	     for (idx in 1:nrow(z.scores)) {
	        names_heatmap_dep = dep.id[idx]
     	        new_names[idx] = sprintf("%i: %s",  					  idx,names_heatmap_dep)
	     }
	  }
	  if(file.exists(sprintf("%s/data/dep.symb.data",run.name))){  
	     symb_all = dget(sprintf("%s/data/dep.symb.data",run.name))
	     symb = symb_all[indices.dep]
     	     for (idx in 1:nrow(z.scores)) {
	        names_heatmap_dep = dep.id[idx]
     	        new_names[idx] = sprintf("%i: %s (%s)", 					  idx,names_heatmap_dep, symb[idx])
	     }
     	  }
        new_names[p.values > significance] <- ""
        rownames(z.scores) = new_names

      }

      if (show.names.indep) {
	  if(!file.exists(sprintf("%s/data/indep.id.data",run.name)))
   	  {
		 msg = paste("please first run the read.data() function to read in the ID's of the data necessary for this function")
    	 	 stop(msg)
   	  }
	  indep.id_all <- dget(sprintf("%s/data/indep.id.data",run.name))
	  indep.id <- indep.id_all[indices.indep]
        new_names = rep("", ncol(z.scores))

	  if(!file.exists(sprintf("%s/data/indep.symb.data",run.name)))
	  {
	     for (idx in 1:ncol(z.scores)) {
	        names_heatmap_indep = indep.id[idx]
     	        new_names[idx] = sprintf("%i: %s",  					  idx,names_heatmap_indep)
	     }
	  }
	  if(file.exists(sprintf("%s/data/indep.symb.data",run.name))){  
	     symb_all = dget(sprintf("%s/data/indep.symb.data",run.name))
	     symb = symb_all[indices.indep]
     	     for (idx in 1:ncol(z.scores)) {
	        names_heatmap_indep = indep.id[idx]
     	        new_names[idx] = sprintf("%i: %s (%s)", 					  idx,names_heatmap_indep, symb[idx])
	     }
     	  }
        new_names[zscores.vec == 1] <- ""
        colnames(z.scores) = new_names

      }
      # Open the graphics device.
    if(pdf==TRUE){
      if ( dev.cur() == 1 ) {
        device.opened = TRUE
        filename = sprintf("%s/heatmap_zscores/globaltest_heatmap_%s_%s_pvals_%s.pdf",run.name,input.region.name,adjust.method,plot.method)
	  A4(filename)

      } else {
        device.opened = FALSE

      }
    }
    else
	device.opened = FALSE


      par(oma=c(1,1,1, 1))

      # Set defaults.
      def.title  = input.region.name
      def.colors = maPalette(high="green", mid="black", low="red")
      def.legend = pval.legend
      def.fill   = vec.col
      if(length(scale) == 1) scale = c(-zmax, zmax)
      # Create the first plot: no clustering whatsoever.

       heatmap.with.plot(dep.data, input.region = input.region.name, plot.method = plot.method, Normal.data = Normal.data, subtype, acgh.heatmap.scale = acgh.heatmap.scale, run.name = run.name, windowsize, lambda=lambda, x=z.scores, Rowv=NA, Colv=NA, RowSideColors=vec.col.pvals,                                         
              ColSideColors=vec.col.top, scale="none", col=def.colors,                                                   
               main=def.title, mar=c(10, 10), zlim=scale, windowsize = windowsize,...)

      image.plot(legend.only=TRUE, zlim=scale,
                 col=maPalette(high="green", low="red", mid="black"))
    par(mar = c(0, 0, 0, 0))
           legend("topright",legend=def.legend, cex=0.5, fill=def.fill)

      # All done. Close the graphics device.

      if (device.opened)
        dev.off()


  if (device.opened) {

       writeln()
       writeln("**************************************************")
       writeln("The results are stored in the following file:")
       writeln(filename)
       writeln("**************************************************")
       writeln()

    }
}
  }
}

