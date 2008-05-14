\name{integrated.analysis}
\alias{integrated.analysis}
\title{Integrated analysis of expression and copy number microarray data}
\description{Runs the Integrated Analysis to test for associations between DNA copy number measurements and gene expression measurements on the same set of samples.}

\usage{integrated.analysis(samples, input.regions="all chrs", adjust=FALSE, zscores=FALSE, 
  method=c("auto", "asymptotic", "permutations", "gamma"), run.name=NULL)
}

\arguments{
 \item{samples}{\code{\link{vector}} with either the names of the columns in the dependent and independent     data corresponding to the samples, or a numerical vector containing the column numbers to include in the analysis, e.g. 5:10 means columns 5 till 10. Make sure that both datasets have the same number of samples with the same column names!}
      
  \item{input.regions}{\code{\link{vector}} indicating the regions to be analyzed. Can be defined in four ways:
  \code{1) predefined input region: } insert a predefined input region, choices are: \code{"all chrs"}, \code{"all chrs auto"}, \code{"all arms"}, \code{"all arms auto"} In the predefined regions \code{"all arms"} and \code{"all arms auto"} the arms 13p, 14p, 15p, 21p and 22p are left out, because in most studies there are no or few probes in these regions. To include them, just make your own \link{vector} of arms. \code{2) whole chromosome(s): }insert a single chromosome or a list of chromosomes as a \code{\link{vector}:} \code{c(1, 2, 3)}. \code{3) chromosome arms: } insert a single chromosome arm  or a list of chromosome arms like \code{c("1q", "2p", "2q")}.\code{4) subregions of a chromosome: } insert a chromosome number followed by the start and end position like \code{c("chr1\_1-1000000")} These regions can also be combined, e.g. \code{c("chr1_1-1000000","2q", 3)}. See \code{details} for more information.}
    
  \item{adjust}{Confounders for which the integration-test must be adjusted, such as tumor type, location, gender, etc. Either a formula with a factor of a vector with names with the same length as samples or FALSE. A formula with a vector can e.g. be: \code{Y~factor(subtype)} where subtype is a vector with the same length as samples with names like: subtype = c("tumor","tumor", "normal","normal", etc...) See \code{?globaltest} for more information.}

	\item{zscores}{\code{Boolean}, indicates whether the z-scores are calculated (takes longer time to run). If z-scores=FALSE, only p-values are calculated.}

  \item{method}{The method for calculation of the p-values. Use \code{method = "asymptotic"} for the full asymptotic distribution of the test statistic, \code{method = "gamma"} for the gamma (= scaled chi-squared) approximation to that distribution and \code{method = "permutations"} for a permutation p-value. The recommended default: \code{method = "auto"} chooses the permutations method if the number of possible permutations does not exceed 10,000 and the asymptotic otherwise. See \code{?globaltest} for more information.}

  \item{run.name}{Name of the analysis. The results will be stored in a folder with this name in the current working directory (use \code{getwd()} to print the current working directory). If the \code{run.name = NULL}, the default folder \code{"analysis\_results"} will be generated.} 
}
\details{
 	The Integrated Analysis is a regression of the independent data
  on the dependent features. In most cases, the dependent data will be the copy 
  number measurements from array-CGH and the independent data the expression array 
  values. The regression itself is done using the \code{globaltest}, which means 
  that the genes in a region (e.g. a chromosome arm) are tested as a gene set. The 
  individual associations between each copy number probe and each expression probe 
  are calculated as z-scores (standardized influences, see \code{?globaltest}).
  
  This function splits the datasets into separate sets for each region (as 
  specified by the \code{input.regions}) and runs the analysis for each region
  separately.
  
  When running the Integrated Analysis for a predefined input region, like \code{"all arms"} 
  and \code{"all chrs"}, output can be obtained for all input regions, as well as 
  subsets of it. But note that the genomic unit must be the same: if \code{integrated.analysis} 
  was run using chromosomes as units, any of the functions and plots must also use chromosomes 
  as units, and not chromosome arms. Similarly, if \code{integrated analysis} was run using 
  chromosome arms as units, these units must also be used to produce plots and outputs. 
  For example if the \code{input.regions = "all arms"} was used, p-value plots 
  (see \link{sim.plot.pvals.on.region} can be produced by inserting the \code{input.regions = "all arms"}, 
  but also for instance \code{"1p"} or \code{"20q"}. However, to produce a plot of the whole 
  chromosome, for example chromosome 1, the integrated should be re-run with \code{input.region=1}.
  The same goes for "all chrs": p-value plots etc. can be produced for chromosome 1,2 and so on... 
  but to produce plots for an arm, the \code{integrated.analysis} should be re-run for that region. 
  This also goes for subregions of the chromosome like "chr1\_1-1000000". 

}
\value{
 	No values are returned. Instead, the results of the analysis are stored
  in the subdirectories of the directory specified in \code{run.name}. E.g. the z-score matrices 
  are saved in subfolder intermediate.data. The following functions can be used to
  visualize the data:
  
  \item{1) }{\code{\link{sim.plot.zscore.heatmap}} (pdf, only possible when \code{zscores=TRUE})}
  \item{2) }{\code{\link{sim.plot.pvals.on.region}} (pdf)}
  \item{3) }{\code{\link{sim.plot.pvals.on.genome}} (pdf)}

  \item{}{Other functions can be used to tabulate the results:}
  \item{1) }{\code{\link{tabulate.pvals}} (data.frame)}
  \item{2) }{\code{\link{tabulate.top.dep.features}} (txt)}
  \item{3) }{\code{\link{tabulate.top.indep.features}} (txt, only possible when zscores=TRUE}
}

\references{\item{1}{Goeman JJ, van de Geer SA, de Kort F, van Houwelingen HC. 
A global test for groups of genes: testing association with a clinical outcome. Bioinformatics. 2004; 20:93-109}}


\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\code{\link{SIM}}, \code{\link{assemble.data}}, \code{\link{sim.plot.zscore.heatmap}}, 
\code{\link{sim.plot.pvals.on.region}}, \code{\link{sim.plot.pvals.on.genome}}, \code{\link{tabulate.pvals}}, 
\code{\link{tabulate.top.dep.features}}, \code{\link{tabulate.top.indep.features}}, 
\code{\link{impute.nas.by.surrounding}}, \code{\link{sim.update.chrom.table}},\code{\link[globaltest]{globaltest}}
}

\examples{
#load the datasets and the samples to run the integrated analysis
data(expr.data)
data(acgh.data)
data(samples) 
         
#assemble the data
assemble.data(dep.data=acgh.data, indep.data=expr.data, ann.dep=colnames(acgh.data)[1:4], ann.indep= colnames(expr.data)[1:4], dep.id="ID", dep.chr="CHROMOSOME", dep.pos="STARTPOS", dep.symb="Symbol", indep.id="ID", indep.chr="CHROMOSOME", indep.pos="STARTPOS", indep.symb="Symbol", overwrite=TRUE, run.name="chr8")

#run the integrated analysis
integrated.analysis(samples=samples, input.regions=c(8), adjust=FALSE, zscores=TRUE, method="auto", run.name ="chr8")

}
\keyword{misc}

