`adj.bhfdr` <-
function(pvals, chrn=0) {
  #if the length of pvalue-vector is greater than 0, the pvalues will be adjusted
  if (length(pvals) > 1) {
  
    # Compute the adjusted p-values using Benjamini & Hochberg correction.
    adj.pvals = mt.rawp2adjp(pvals, "BH")
    
    # Order the adjusted p-values on their *original* index.
    adj.pvals = adj.pvals$adjp[order(adj.pvals$index), 2]

  } else {
    adj.pvals = pvals
    
  }  
  # Return the adjusted p-values.
  return(adj.pvals)
}

