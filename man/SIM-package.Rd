\name{SIM-package}
\alias{SIM-package}
\alias{SIM}
\docType{package}
\title{Statistical Integration of Microarrays}
\description{SIM is a statistical model to identify copy number changes that affect the expression 
of genes within the same chromosomal region. Copy number is considered as the dependent 
variable and expression as the independent variable. Copy number alterations may span 
many expression probes and affect them in a possibly subtle but consistent way.
Therefore, we test whether copy number is associated with a set of expression levels 
within a chromosome arm (or mimimal common region) in a random-effect model. Association scores for individual 
expression levels (z-scores) are also calculated. For more information on the 
random-effect model, see \code{?globaltest}.

Each sample should be profiled both on a copy number and on an expression array. The array 
platforms used for DNA and RNA analysis may be different as long as the probes have mapped 
to the genome. RESOURCERER can be used to search chromosome and basepair location for 
expression microarray probes \code{(http://compbio.dfci.harvard.edu/tgi/cgi-bin/magic/r1.pl)}. 
See \code{RESOURCERER.annotation.to.ID} on how to insert this information as annotation columns.
Alternatively, the chromosome, basepair locations and gene symbol can be extracted from 
AnnotationData packages available in Bioconductor or generated using the AnnBuilder package.

When copy number data is run as dependent variable, we use \code{method.adjust="BY"} for 
multiple testing correction. This method accounts for dependence between measurements and 
is more conservative than "BH". For details on the multiple testing correction methods see 
?p.adjust. We have experienced that a rather low stringency cut-off on the BY-values of 
20\% allows the detection of associations for data with a low number of samples or a low 
frequency of abberations. False positives are rarely observed.

Make sure that the array probes are mapped to the same builds of the genome, and that the 
\code{chrom.table} used by the \code{\link{integrated.analysis}} is from the same build as well.
See \code{\link{sim.update.chrom.table}}.
}
\details{
\tabular{ll}{
Package: \tab SIM\cr
Type: \tab Package\cr
Version: \tab 1.9.0\cr
Date: \tab 2008-02-06\cr
License: \tab Open\cr
}
}
 
\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\references{R.X. de Menezes, M. Boetzer, M. Sieswerda, G.J.B. van Ommen, J.M. Boer
Integrated Statistical analysis to identify associations between DNA copy number and gene expression in microarray data. Submitted.}

\keyword{package}

\seealso{
\code{\link{assemble.data}}, \code{\link{integrated.analysis}}, \code{\link{sim.plot.zscore.heatmap}}, 
\code{\link{sim.plot.pvals.on.region}}, \code{\link{sim.plot.pvals.on.genome}}, \code{\link{tabulate.pvals}}, 
\code{\link{tabulate.top.dep.features}}, \code{\link{tabulate.top.indep.features}}, 
\code{\link{impute.nas.by.surrounding}}, \code{\link{sim.update.chrom.table}}
}

\examples{
#load the datasets and the samples to run the integrated analysis
data(expr.data)
data(acgh.data)
data(samples) 
         
#assemble the data
assemble.data(dep.data = acgh.data, indep.data = expr.data,ann.dep = colnames(acgh.data)[1:4], ann.indep = colnames(expr.data)[1:4], dep.id="ID", dep.chr = "CHROMOSOME",dep.pos = "STARTPOS",dep.symb="Symbol",  indep.id="ID",indep.chr = "CHROMOSOME", indep.pos = "STARTPOS", indep.symb="Symbol", overwrite = TRUE,run.name = "chr8")

#run the integrated analysis
integrated.analysis(samples = samples, input.regions = 8, adjust=FALSE, zscores=TRUE, method = "auto", run.name = "chr8")

# use functions to plot the results of the integrated analysis

#plot the p-values along the genome
sim.plot.pvals.on.genome(input.regions = 8,adjust.method = "BY",pdf = FALSE, run.name = "chr8")

#plot the p-values along the regions
sim.plot.pvals.on.region(input.regions = 8, adjust.method="BY", run.name = "chr8")

#plot the z-scores in an association heatmap
sim.plot.zscore.heatmap(input.regions = 8, significance=0.2, z.threshold=3, show.names.dep=TRUE,show.names.indep=TRUE, adjust.method = c("BY"), scale="auto", plot.method = "smooth", pdf = FALSE, run.name = "chr8")

#tabulate the p-values per region (prints to screen)
tabulate.pvals(input.regions = 8,adjust.method="BY", bins=c(0.001,0.005,0.01,0.025,0.05,0.075,0.10,0.20,1.0), significance.idx=8, order.by="\%", decreasing=TRUE, run.name = "chr8") 

#get the top dependent features sorted by p-value
tabulate.top.dep.features(input.regions = 8, adjust.method="BY",run.name = "chr8")

#get the top independent features sorted by mean z-score
tabulate.top.indep.features(input.regions = 8,adjust.method="BY", significance=0.2, sort.order='positive', run.name = "chr8")
}


