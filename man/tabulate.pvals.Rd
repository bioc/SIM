\name{tabulate.pvals}
\alias{tabulate.pvals}
\title{Sums significant p-values for the analyzed regions}
\description{Generates a \code{\link{data.frame}} with the signicance of p-values in the analyzed regions, dividing them into bins.
}

\usage{tabulate.pvals(input.regions = "all chrs", adjust.method = "BY", bins = c(0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 1), significance.idx = 8, order.by = "\%", decreasing = TRUE, run.name = NULL)
}
\arguments{
  \item{input.regions}{\code{\link{vector}} with analyzed regions for which to produce the table. 
  Can be defined in four ways:

  \code{1) predefined input region: } insert a predefined input region, choices are:
  \code{"all chrs"}, \code{"all chrs auto"}, \code{"all arms"}, \code{"all arms auto"}
  In the predefined regions \code{"all arms"} and \code{"all arms auto"} the arms 13p, 
  14p, 15p, 21p and 22p are left out, because in most studies there are no or few probes 
  in these regions. To include them, just make your own \code{\link{vector}} of arms.
  \code{2) whole chromosome(s): }insert a single chromosome or a list of chromosomes as a \code{\link{vector}:
}
  \code{c(1, 2, 3)}.
  \code{3) chromosome arms: } insert a single chromosome arm  or a list of chromosome arms like 
  \code{c("1q", "2p", "2q")}.
  \code{4) subregions of a chromosome: } insert a chromosome number followed by the start and end position lik
e \code{c("chr1_1-1000000")}
  These regions can also be combined, e.g. \code{c("chr1_1-1000000","2q", 3)}.
  See \code{details} for more information.}  
  
  \item{adjust.method}{Method used to adjust the p-values for multiple
   testing. Either \code{"BY"} (recommended when copy number is used as dependent data), 
   \code{"BH"} or \code{"raw"}. See \code{\link{SIM}} for more information about adjusting p-values.}
     
  \item{bins}{\code{\link{vector}} of significance thresholds. This function will calculate 
  the number of features having a p-value lower than the bin.}

  \item{significance.idx}{Index of "bins" to use when computing the percentage
   of significant p-values. Defaults to 8 (i.e. the first entry in "bins"), in this case 0.20.}
  
  \item{order.by}{Column used for sorting the table. Defaults to "\%" (i.e. the percentage of significant p-va
lues).}
  
  \item{decreasing}{Direction used for sorting.
   Defaults to TRUE (i.e. highest values on top).}

  \item{run.name}{Name of the analysis. The results will be 
   stored in a folder with this name in the current working directory 
   (use \code{getwd()} to print the current working directory). 
   If the \code{run.name = NULL}, the default folder \code{"analysis_results"} will be generated.}
}


\value{
 Returns a \code{data.frame}. Each row corresponds to a chromosome and has
  as many entries as entries in bins, plus 1. Each entry contains the 
  number of p-values that is smaller or equal to the corresponding entry
  in bins. 
  
  The last entry holds the percentage of p-values that is smaller than or
  equal to the bin identified by \code{significance.idx}.
}

\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\code{\link{SIM}}, \code{\link{assemble.data}}, \code{\link{integrated.analysis}}, 
\code{\link{sim.plot.zscore.heatmap}}, \code{\link{sim.plot.pvals.on.region}}, 
\code{\link{sim.plot.pvals.on.genome}}, \code{\link{tabulate.top.dep.features}}, 
\code{\link{tabulate.top.indep.features}}, \code{\link{impute.nas.by.surrounding}}, 
\code{\link{sim.update.chrom.table}}
}


\examples{
#load the datasets and the samples to run the integrated analysis
data(expr.data)
data(acgh.data)
data(samples) 
         
#assemble the data
assemble.data(dep.data = acgh.data, indep.data = expr.data, ann.dep = 1:4,ann.indep = 1:4, dep.id="ID", dep.chr = "CHROMOSOME",dep.pos = "STARTPOS",dep.symb="Symbol", indep.id="ID",indep.chr = "CHROMOSOME", indep.pos = "STARTPOS", indep.symb="Symbol", overwrite = TRUE,run.name = "chr8")

#run the integrated analysis
integrated.analysis(samples = samples, input.regions = 8, adjust = FALSE, zscores=TRUE, method = "auto",run.name = "chr8")

#tabulate the p-values per region
tabulate.pvals(input.regions = 8,adjust.method="BY", bins=c(0.001,0.005,0.01,0.025,0.05,0.075,0.10,0.20,1.0), significance.idx=8, order.by="\%", decreasing=TRUE, run.name = "chr8") 

}
\keyword{misc}

