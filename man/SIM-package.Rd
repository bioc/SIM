\name{SIM-package}

\alias{SIM-package}

\alias{SIM}

\docType{package}

\title{Statistical Integration of Microarrays}

\description{
	SIM is a statistical model to identify associations between two genomic datasets. Where one is assigned as
dependent variable and the other as independent e.g. copy number measurements on several samples versus 
expression measurements on the same samples. 
A region of interest can be chosen to run the integrated analysis on either the same region for both dependent and 
independent datasets or different regions. For each dependent feature a P-value measures the association with the
independent data, the contribution of each independent feature is given as Z-scores. The integrated analysis 
is based on the random-effect model for gene-sets as implemented in \link[globaltest:v4_globaltest]{globaltest}.

maybe something about annotation?

By default we use \code{method.adjust = "BY"} (Benjamini-Yekutieli) for multiple testing correction. 
This method accounts for dependence between measurements and is more conservative than "BH" (Benjamini-Hochberg). 
For details on the multiple testing correction methods see \link[stats:p.adjust]{p.adjust}. 
We have experienced that a rather low stringency cut-off on the BY-values of 
20\% allows the detection of associations for data with a low number of samples or a low 
frequency of abberations. False positives are rarely observed.

Make sure that the array probes are mapped to the same builds of the genome, and that the 
\link{chrom.table} used by the \link{integrated.analysis} is from the same build as well.
See \link{sim.update.chrom.table}.
}

\details{
\tabular{ll}{
Package: \tab SIM\cr
Type: \tab Package\cr
Version: \tab 1.99.0\cr
Date: \tab 2009-08-13\cr
License: \tab Open\cr}
}
 
\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\references{
Menezes RX, Boetzer M, Sieswerda M, van Ommen GJ, Boer JM (2009). 
Integrated analysis of DNA copy number and gene expression microarray data using gene sets.
\emph{BMC Bioinformatics}, \bold{10}, 203-.

Goeman JJ, van de Geer SA, de Kort F, van Houwelingen HC (2004). 
A global test for groups of genes: testing association with a clinical outcome. 
\emph{Bioinformatics}, \bold{20}, 93-109.
}

\keyword{package}

\seealso{
\link{assemble.data}, 
\link{integrated.analysis}, 
\link{sim.plot.zscore.heatmap}, 
\link{sim.plot.pvals.on.region}, 
\link{sim.plot.pvals.on.genome}, 
\link{tabulate.pvals}, 
\link{tabulate.top.dep.features}, 
\link{tabulate.top.indep.features}, 
\link{impute.nas.by.surrounding}, 
\link{sim.update.chrom.table}, 
\link{sim.plot.overlapping.indep.dep.features},
\link{getoverlappingregions}
}
\examples{
#load the datasets and the samples to run the integrated analysis
data(expr.data)
data(acgh.data)
data(samples) 
         
#assemble the data
assemble.data(dep.data = acgh.data, 
              indep.data = expr.data,
              dep.ann = colnames(acgh.data)[1:4], 
              indep.ann = colnames(expr.data)[1:4], 
              dep.id="ID", 
              dep.chr = "CHROMOSOME",
              dep.pos = "STARTPOS",
              dep.symb="Symbol",  
              indep.id="ID",
              indep.chr = "CHROMOSOME", 
              indep.pos = "STARTPOS", 
              indep.symb="Symbol", 
              overwrite = TRUE,
              run.name = "chr8q")

#run the integrated analysis
integrated.analysis(samples = samples, 
                    input.regions ="8q", 
                    zscores=TRUE, 
                    run.name = "chr8q")

# use functions to plot the results of the integrated analysis

#plot the P-values along the genome
sim.plot.pvals.on.genome(input.regions = "8q", 
                         significance = c(0.2, 0.05), 
                         adjust.method = "BY", 
                         pdf = FALSE, 
                         run.name = "chr8q")

#plot the P-values along the regions
sim.plot.pvals.on.region(input.regions = "8q", 
						 adjust.method="BY", 
						 run.name = "chr8q")

#plot the z-scores in an association heatmap
#plot the zscores in a heatmap
sim.plot.zscore.heatmap(input.regions = "8q", 
                        method="full", 
                        significance=0.2,                        
                        z.threshold=3, 
                        show.names.indep=TRUE, 
                        show.names.dep=TRUE, 
                        adjust.method = "BY",  
                        add.plot = "smooth", 
                        smooth.lambda = 2,
                        pdf = FALSE, 
                        run.name = "chr8q")

sim.plot.zscore.heatmap(input.regions = "8q", 
                        method="full", 
                        significance = 0.05,                        
                        z.threshold = 1, 
                        show.names.indep=TRUE, 
                        show.names.dep=FALSE, 
                        adjust.method = "BY",  
                        add.plot = "heatmap", 
                        smooth.lambda = 2,
                        pdf = FALSE, 
                        run.name = "chr8q")
                        
sim.plot.zscore.heatmap(input.regions = "8q", 
                        method="full", 
                        significance = 0.05,                        
                        z.threshold = 1, 
                        show.names.indep=TRUE, 
                        show.names.dep=TRUE, 
                        adjust.method = "BY",  
                        add.plot = "none",                                                                      
                        pdf = FALSE, 
                        run.name = "chr8q")

#tabulate the P-values per region (prints to screen)
tabulate.pvals(input.regions = "8q", 
               adjust.method="BY", 
               bins=c(0.001,0.005,0.01,0.025,0.05,0.075,0.10,0.20,1.0),               
               run.name = "chr8q") 
               
table.dep <- tabulate.top.dep.features(input.regions="8q", 
		                  adjust.method="BY", 
						  method="full",
						  significance=0.05,
						  run.name="chr8q")
head(table.dep[["8q"]])

table.indep <- tabulate.top.indep.features(input.regions="8q",
		                                  adjust.method="BY",
										  method="full",
										  significance= 0.05,										  
										  z.threshold=c(-1, 1),
										  run.name="chr8q")
head(table.indep[["8q"]])	

sim.plot.overlapping.indep.dep.features(input.regions="8q", 
		                                adjust.method="BY", 
										significance=0.1, 
										z.threshold= c(-1,1),                                        
										log=TRUE,	
										summarize="consecutive",										
										pdf=FALSE, 
										method="full",
										run.name="chr8q") 									  
						  
}

