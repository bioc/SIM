\name{sim.plot.pvals.on.genome}
\alias{sim.plot.pvals.on.genome}
\title{Plot the p-values in whole genome overview}
\description{Generates a plot of the analyzed dependent data probe positions and their significance on all chromosomes.}

\usage{
sim.plot.pvals.on.genome(input.regions="all chrs", adjust.method=c("BY", "BH", "raw"), pdf = TRUE, run.name=NULL,...)
}

\arguments{
 \item{input.regions}{\code{\link{vector}} with analyzed regions for which to produce the graphs. Can be defined in four ways: \code{1) predefined input region: } insert a predefined input region, choices are: \code{"all chrs"}, \code{"all chrs auto"}, \code{"all arms"}, \code{"all arms auto"} In the predefined regions \code{"all arms"} and \code{"all arms auto"} the arms 13p, 14p, 15p, 21p and 22p are left out, because in most studies there are no or few probes in these regions. To include them, just make your own \link{vector} of arms. \code{2) whole chromosome(s): }insert a single chromosome or a list of chromosomes as a \code{\link{vector}:} \code{c(1, 2, 3)}. \code{3) chromosome arms: } insert a single chromosome arm  or a list of chromosome arms like \code{c("1q", "2p", "2q")}. \code{4) subregions of a chromosome: } insert a chromosome number followed by the start and end position like \code{c("chr1_1-1000000")} These regions can also be combined, e.g. \code{c("chr1_1-1000000","2q", 3)}. For more information see the \code{details} section of \link{integrated.analysis}.}  
  
  \item{adjust.method}{Method used to adjust the p-values for multiple testing. Either \code{"BY"} (recommended when copy number is used as dependent data), \code{"BH"} or \code{"raw"}. Defaults to "BY". See \link{SIM} for more information about adjusting p-values.}
\item{pdf}{Logical. Indicate whether to generate a pdf of the plot in the run.name directory or plot on screen.} 
	\item{run.name}{Name of the analysis. The results will be stored in a folder with this name in the current working directory (use \code{getwd()} to print the current working directory). If the \code{run.name = NULL}, the default folder \code{"analysis_results"} will be generated.}
\item{...}{Arguments to be passed to methods, such as graphical parameters (see \link{par}).}
}

\details{
The orange triangles indicate the start and end of the analyzed regions. The purple dot indicates the cent
romere.}

\value{No values are returned. The results are stored in the folder \code{run.name} as pdf.}
\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{\code{\link{SIM}}, \code{\link{assemble.data}}, \code{\link{integrated.analysis}}, \code{\link{sim.plot.zscore.heatmap}}, \code{\link{sim.plot.pvals.on.region}}, \code{\link{tabulate.pvals}}, \code{\link{tabulate.top.dep.features}}, \code{\link{tabulate.top.indep.features}}, \code{\link{impute.nas.by.surrounding}}, \code{\link{sim.update.chrom.table}}
}

\examples{
#load the datasets and the samples to run the integrated analysis
data(expr.data)
data(acgh.data)
data(samples) 
         
#assemble the data
assemble.data(dep.data=acgh.data, indep.data=expr.data, ann.dep=colnames(acgh.data)[1:4], ann.indep=colnames(expr.data)[1:4], dep.id="ID", dep.chr="CHROMOSOME", dep.pos="STARTPOS", dep.symb="Symbol", indep.id="ID", indep.chr="CHROMOSOME", indep.pos="STARTPOS", indep.symb="Symbol", overwrite=TRUE, run.name="chr8")

#run the integrated analysis
integrated.analysis(samples=samples, input.regions=c(8), adjust=FALSE, zscores=TRUE, method="auto", run.name ="chr8")

#plot the p-values along the genome
sim.plot.pvals.on.genome(input.regions=c(8), adjust.method="BY", pdf=FALSE, run.name="chr8")
}
\keyword{misc}

