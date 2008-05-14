`sim.plot.adj.pvals` <-
function(input.region.name, adjust.method=c("BY", "BH","raw"), run.name = NULL) {
  # Open device.
    options(warn=-1)
    raw.pvals  = try(get.pvals(input.region.name,run.name = run.name), silent=T)
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

      ngenes = length(raw.pvals)

      if (ngenes > 1) {
        par(mfrow=c(1,3))
      
        # Histogram.
        hist(raw.pvals, main="Raw p-values", xlab="", col="blue")
        

        # Plot of raw p-values.
        subtitle = sprintf("Input region: %s", input.region.name)

        m = "Sorted raw p-values"
        m = sprintf("Sorted raw p-values for %s", input.region.name)
        plot(1:ngenes, sort(raw.pvals), main=m, xlab="Genes",
             ylab="", pch=20, col="blue", ylim=c(0,1), sub=subtitle)

        segments(0, 0, ngenes, 1, lty="dashed")

        # Plot of adjusted p-values.
        m = sprintf("%s-corrected p-values", adjust.method)

        plot(1:ngenes, sort(p.values), main=m, xlab="Genes", 
             ylab="Sorted adjusted p-values", pch=20, col="blue", ylim=c(0,1))

        segments(0, 0.2, ngenes, 0.2, lty="dashed")
      }

    } 
}

