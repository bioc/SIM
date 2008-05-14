\name{tabulate.top.indep.features}
\alias{tabulate.top.indep.features}
\title{Lists the mean z-scores for the independent features}
\description{Lists the mean z-scores for independent features in the analyzed regions, calculated 
across the significant dependent features. Gives insight in the expression levels most 
strongly associated with copy number changes.
}

\usage{
tabulate.top.indep.features(input.regions = "all chrs", adjust.method = c("BY", "BH", "raw"), significance = 0.2, sort.order = "positive", run.name = NULL)
}

\arguments{

  \item{input.regions}{\code{\link{vector}} indicating the regions to be analyzed. Can be defined in four 
ways:
  \code{1) predefined input region: } insert a predefined input region, choices are:
  \code{"all chrs"}, \code{"all chrs auto"}, \code{"all arms"}, \code{"all arms auto"}
  In the predefined regions \code{"all arms"} and \code{"all arms auto"} the arms 13p, 
  14p, 15p, 21p and 22p are left out, because in most studies there are no or few probes 
  in these regions. To include them, just make your own \code{\link{vector}} of arms.
  \code{2) whole chromosome(s): }insert a single chromosome or a list of chromosomes as a \code{\link{vector}}
  \code{c(1, 2, 3)}.
  \code{3) chromosome arms: } insert a single chromosome arm  or a list of chromosome arms like 
  \code{c("1q", "2p", "2q")}.
  \code{4) subregions of a chromosome: } insert a chromosome number followed by the start and end position
 like \code{c("chr1_1-1000000")}
  These regions can also be combined, e.g. \code{c("chr1_1-1000000","2q", 3)}.
  See \code{details} for more information.}

  \item{adjust.method}{Method used to adjust the p-values for multiple
  testing. Either \code{"BY"} (recommended when copy number is used as dependent data), 
  \code{"BH"} or \code{"raw"}. Defaults to "BY". See \code{\link{SIM}} for more information about adjustin
g p-values.}
   
  \item{significance}{threshold used to select the significant dependent
   features. Only the z-scores with these features are used to calculate the mean z-scores 
   across the independent features.}

  \item{sort.order}{ Indicates how the z-scores are sorted, either \code{"positive"} or \code{"negative"}.}

  \item{run.name}{Name of the analysis. The results will be stored in a folder with this 
  name in the current working directory (use \code{getwd()} to print the current working 
  directory). If the \code{run.name = NULL}, the default folder \code{"analysis_results"} 
  will be generated.} 
}

\details{
\code{tabulate.top.indep.features} can only be run after \code{\link{integrated.analysis}} 
with \code{zscores=T}.

Output is a .txt file containing a table with the mean z-scores of all independent features 
per analyzed region. It includes the \code{ann.indep} columns that were read in the 
\code{\link{assemble.data}} function.

Depending on the argument "adjust.method", the p-values are first corrected for multiple testing. Next, th
e z-scores are filtered to include only those entries that correspond to significant (p-value < "significa
nce") dependent features to calculate the mean z-scores.  
}
\value{No values are returned. The results are stored in a subdirectory of \code{run.name} as pdf.}

\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\code{\link{SIM}}, \code{\link{assemble.data}}, \code{\link{integrated.analysis}}, 
\code{\link{sim.plot.zscore.heatmap}}, \code{\link{sim.plot.pvals.on.region}}, 
\code{\link{sim.plot.pvals.on.genome}}, \code{\link{tabulate.pvals}}, 
\code{\link{tabulate.top.dep.features}}, \code{\link{impute.nas.by.surrounding}}, 
\code{\link{sim.update.chrom.table}}
}



\examples{
#load the datasets and the samples to run the integrated analysis
data(expr.data)
data(acgh.data)
data(samples) 
         
#assemble the data
assemble.data(dep.data = acgh.data, indep.data = expr.data, ann.dep = colnames(acgh.data)[1:4], ann.indep = colnames(expr.data)[1:4], dep.id="ID", dep.chr = "CHROMOSOME",dep.pos = "STARTPOS",dep.symb="Symbol",indep.id="ID",indep.chr = "CHROMOSOME", indep.pos = "STARTPOS",indep.symb="Symbol", overwrite = TRUE,run.name = "chr8")

#run the integrated analysis
integrated.analysis(samples = samples, input.regions = 8, adjust=FALSE, zscores=TRUE, method = "auto", run.name="chr8")

#get the highest associated independent features
tabulate.top.indep.features(input.regions = 8, adjust.method="BY", significance=0.2,sort.order='positive', run.name = "chr8")
}
\keyword{misc}
