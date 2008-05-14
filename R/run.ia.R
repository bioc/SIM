`run.ia` <-
function(indep.chr.matrix, dep.chr.matrix, subtype, zscores, method =  c("auto", "asymptotic", "permutations", "gamma")) {
  
  #
  # Create matrices where results of the globaltest will be stored in.
  #
  # Create a matrix to hold, for each dependent feature, the processed z-scores that
  # indicate how each independent features influences the dependent feature's p-value.
  # Rows correspond to dependent features.
  gene.proc.zscores = array(0, dim = c(nrow(dep.chr.matrix), nrow(indep.chr.matrix)))
  
  # Create a matrix to hold colors. "1" indicates a positive correlation
  # with a feature and "2' indicates a negative correlation.
  # Rows correspond to dependent features.
  gene.proc.zscores.col = array(0, dim = c(nrow(dep.chr.matrix), nrow(indep.chr.matrix)))
  

  
  # List to hold one p-value for each feature.
  p.values = rep(0, nrow(dep.chr.matrix))
  
  # Create variables used for estimating the time to completion.
  time.left = Inf
  time.single.iteration = Inf
    
  #
  # Run the globaltest for every aCGH feature.
  #
  
  nrows = nrow(dep.chr.matrix)
  time.start = proc.time()[3]
  
  for(idx in 1:nrows) {
    
    # Determine the time left.
    time.left = time.single.iteration * (nrows - idx)
    
    # Store the current time ...
    time.start.iteration = proc.time()[3]
    
  
    # Print progress every 10 features.
    if (idx %% 10 == 0) {
      time.left.minutes = floor(time.left / 60)
      time.left.seconds = time.left - (time.left.minutes*60)
      
      writeln("  feature %i/%i, time left: %2.0f minutes, %2.0f seconds", 
              idx, nrows, time.left.minutes, time.left.seconds)
    }       
    
    
    # Run the globaltest for one aCGH feature. and calculate the influence of
    # all the genes on it. If "subtype" has more than one level, it is used as
    # a confounder.
    current.feature = dep.chr.matrix[idx, ]
    
    method = match.arg(method, c("auto", "asymptotic", "permutations", "gamma"))

	if(class(subtype) == "formula"){
current.gt = globaltest(indep.chr.matrix, current.feature, 	method=method, adjust=subtype)}
    
     else if(length(subtype) ==1) {
	if(subtype==FALSE){
      current.gt = globaltest(indep.chr.matrix, current.feature, 	method=method)}
    } else{
      if(method == "permutations")
      {
      	 msg = paste("you cannot use permutations when running the analysis as a confounder, please change your method")
   	 stop(msg)
      }
	else if(length(subtype) == ncol(indep.chr.matrix)){
	subtype <- factor(subtype)
      current.gt = globaltest(indep.chr.matrix, current.feature, 	method=method, adjust=Y~subtype)}
else{
current.gt = globaltest(indep.chr.matrix, current.feature, 	method=method, adjust=subtype)}
    }

    
    
    # Get the p-value that indicates whether all of the genes together are
    # significantly associated with the current feature.
    p.values[idx] = p.value(current.gt)
    
    # Determine the z-scores for all genes.
    if (zscores) {

      # Use the method "geneplot" to extract the z-scores.
      current.gp = geneplot(current.gt, drawlabels=FALSE, plot=FALSE)
          
      # Get the z-score for each gene.
      gene.z.scores = z.score(current.gp)
      
      # Set z-scores with an "NA" value and all z-scores < 0 to zero.
      gene.z.scores[is.na(gene.z.scores)] = 0
      gene.z.scores[gene.z.scores < 0] = rep(0, length(gene.z.scores[gene.z.scores < 0]))
    
      # Get the color of the plotted bar. "1" indicates a positive correlation
      # with the current feature and "2' indicates a negative correlation.
      my.up = current.gp@res[, 4]
      
      # Create a copy of the z-scores as computed above, and modify it so that
      # entries that have a negative correlation with current feature, will also
      # have a negative z-score.
      processed.z.scores = gene.z.scores
      processed.z.scores[my.up == 2] = -gene.z.scores[my.up == 2]
      
      # Store the processed gene.z.scores.
      gene.proc.zscores[idx, ] = processed.z.scores
      
      # Store the sign of the correlation with the current feature.
      gene.proc.zscores.col[idx, ] = my.up
    }

        
    # Determine the time spent on the current iteration.
    time.single.iteration = proc.time()[3] - time.start.iteration
  }
   
  if (zscores) {
    retval = list(p.values, gene.z.scores,  gene.proc.zscores,  gene.proc.zscores.col)
  }
    
  if (!zscores) {
    retval = list(p.values)
  }
method2 <- function(gt) {
switch (gt@method, 
"No method",
"Gamma approximation",
"Asymptotic distribution",
paste("All", ncol(gt@PermQs), "permutations"),
paste(ncol(gt@PermQs), "random permutations")
)
}
print(paste("method =", method2(current.gt)))
  # Return "retval".
  return(retval)
  }

