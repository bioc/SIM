`sim.plot.pvals.on.region` <-
function(input.regions = c("all chrs"), adjust.method = c("BY", "BH", "raw"), run.name = NULL,...) {
while(names(dev.cur()) == "pdf") dev.off()
  data(chrom.table)
  miss_arms = c("13p", "14p", "15p", "21p", "22p")
  options(warn=-1)
  #if asked for a predefined input region, get the input.regions:
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
  options(warn=1)
  if(is.null(run.name))
  {
	run.name = "analysis_results"
  }

  if (!file.exists(run.name)) {
    msg = paste("The directory", run.name, "does not exist, please insert a correct run.name.")
    stop(msg)
  }  


  # Open the graphics device.
  filename = sprintf("%s/pvalue_plots/pvals_on_region.pdf",run.name)
  pdf(filename, width=11.69, height=8.27)
   for (input.regions.name in input.regions) {
    # Load the p-values from disk.
    options(warn=-1)
    raw.pvals  = try(get.pvals(input.regions.name,run.name = run.name), silent=T)
    options(warn=-1)
    if ( class(raw.pvals) != "try-error" ) {
      
      # Adjust the p-values for multiple testing.
      adjust.method = match.arg(adjust.method, c("BY", "BH", "raw"))

      if (adjust.method == "BY")
        p.values =  adj.byfdr(raw.pvals)
      else if (adjust.method == "BH")
        p.values = adj.bhfdr(raw.pvals)
      else if (adjust.method == "raw")
        p.values = raw.pvals
      
      # Determine the start of the current chromosome.
      start.input.regions = get.start.input.regions(input.regions.name,run.name = run.name) / 10^6
      
      sim.plot.adj.pvals(input.regions.name, adjust.method, run.name = run.name)

      par(mfrow=c(1,1))      
      # Plot the pvalues along the chromosome.
      caption = paste(adjust.method, "-corrected p-values for ", input.regions.name, sep='')
	if(length(start.input.regions) > 1000) cex <- 0.5
	if(length(start.input.regions) < 1000) cex <- 1
      plot(start.input.regions, p.values, col="blue", pch=20, main=caption, ylim=c(0,1), cex = cex,...)
     
      segments(0, 0.10, max(start.input.regions), 0.10, col="gray", lty="dotted")
      segments(0, 0.20, max(start.input.regions), 0.20, col="gray", lty="dotted")
      
    }
     else {
      writeln("Will not generate plot for '%s': it has no features.", input.regions.name)
    }
  }

  dev.off()
  
  writeln()
  writeln("**************************************************")
  writeln("The results are stored in the following file:")
  writeln(filename)
  writeln("**************************************************")
  writeln()
  
}

