`tabulate.pvals` <-
function(input.regions = "all chrs", adjust.method="BY", 
         bins=c(0.001,0.005,0.01,0.025,0.05,0.075,0.10,0.20,1.0), 
         significance.idx=8, order.by="%", decreasing=TRUE, run.name = NULL) 
{
  data(chrom.table)
  miss_arms = c("13p", "14p", "15p", "21p", "22p")
  if(is.null(run.name))
  {
	run.name = "analysis_results"
  }

if (!file.exists(run.name)) {
    msg = paste("The directory", run.name, "does not exist, please insert a correct run.name.")
    stop(msg)
  }

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

  bins = sort(bins)
  retval = data.frame()
  
  for (rowname in input.regions) {
    # Get the p-values for the current chromosome and adjust them for mulitple
    # testing.
    options(warn=-1)
    raw.pvals = try(get.pvals(rowname, run.name), silent=T)
    options(warn=1)

    
    if (class(raw.pvals) == "try-error")
      raw.pvals = c()

    if (adjust.method == "BY")
      p.values =  adj.byfdr(raw.pvals)
    
    else if (adjust.method == "BH")
      p.values = adj.bhfdr(raw.pvals)
    
    else if (adjust.method == "raw")
      p.values = raw.pvals

    
    # Create the bins that. Each bin holds the number of p-values more or 
    # equally significant than its value (as specified by "bins").
    comp.bins = rep(0, length(bins) + 1)
    
    for(i in 1:length(bins))
       comp.bins[i] = sum(p.values<=bins[i]) 
       
    # The last column holds the percentage of significant p-values on the chromosome.
    comp.bins[ length(comp.bins) ] = 100*(comp.bins[significance.idx] / length(p.values))
    
    # Cast comp.bins to an array and add the current set of bins as a row to the
    # matrix that we'll return.
    comp.bins = array(comp.bins, dim=c(1, length(comp.bins)))
    retval = rbind(retval, comp.bins)
  }
  
  retval = cbind(input.regions, retval)
  
  # Cast the matrix to a data frame and set the colnames.
  retval = as.data.frame(retval)
  colnames(retval) = c("input.region", bins, "%")

  # Sort if necessary.
  if (order.by != "") {
    args = list( retval[[order.by]] )
    
    for (i in 2:length(bins))
      args[[i]] = retval[[i]]
    
    args[['decreasing']] = decreasing
    indices = do.call(order, args)
    
    retval = retval[indices, ]
  }

  return(retval)
}

