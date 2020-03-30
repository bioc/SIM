\name{integrated.analysis}
\alias{integrated.analysis}
\title{Integrated analysis of dependent and indepedent microarray data}
\description{Runs the Integrated Analysis to test for associations between dependent and independent microarray data 
             on the same set of samples.}

\usage{integrated.analysis(samples, 
                           input.regions = "all chrs", 
                           input.region.indep = NULL,	 
                           zscores =  FALSE, 
                           method = c("full", "smooth", "window", "overlap"),  
                           dep.end = 1e5, 
                           window = c(1e6, 1e6), 
                           smooth.lambda=2, 
                           adjust = ~1,
                           run.name = "analysis_results", 
                           ...) 
}

\arguments{
	\item{samples}{\code{vector} with either the names of the columns in the dependent and independent data 
	               corresponding to the samples, or a numerical vector containing the column numbers to include in 
	               the analysis, e.g. 5:10 means columns 5 till 10. Make sure that both datasets have the same number 
	               of samples with the same column names!}
    \item{input.regions}{\code{vector} indicating the dependent regions to be analyzed. Can be defined in four ways:
  	      \code{1) predefined input region: } insert a predefined input region, choices are: 
  	      \dQuote{all chrs}, 
  	      \dQuote{all chrs auto}, 
  	      \dQuote{all arms}, 
  	      \dQuote{all arms auto} 
  	      In the predefined regions \dQuote{all arms} and \dQuote{all arms auto} the arms 13p, 14p, 15p, 21p and 22p 
  	      are left out, because in most studies there are no or few probes in these regions. 
  	      To include them, just make your own vector of arms. 
  	      \code{2) whole chromosome(s): } insert a single chromosome or a list of chromosomes as a 
  	      \code{vector:} 
  	      \code{c(1, 2, 3)}. 
  	      \code{3) chromosome arms: } insert a single chromosome arm  or a list of chromosome arms like 
  	      \code{c("1q", "2p", "2q")}.
  	      \code{4) subregions of a chromosome: } insert a chromosome number followed by the start and end position like 
  	      \code{"chr1:1-1000000"} 
  	      These regions can also be combined, e.g. \code{c("chr1:1-1000000","2q", 3)}. 
  	      See \code{details} for more information.}
    \item{input.region.indep}{indicating the independent region which will be analysed in combination of the dependent 
                              region. Only one input region can given using the same format as the dependent input region.}  	       
	\item{zscores}{\code{logical} indicates whether the Z-scores are calculated (takes longer time to run). 
	               If \code{zscores = FALSE}, only P-values are calculated.}
	\item{method}{either one of \dQuote{full}, \dQuote{window}, \dQuote{overlap} or \dQuote{smooth}. This defines how the data is used for the\code{integrated.analysis}.  
		full: the whole dependent data region is taken.
		window: takes the middle of the dependent probe and does the integration on the independent probes that are within 
		        the window given at window-size given by \code{window}.
		overlap: does the integration on the independent probes that are within the start and end of the dependent probes given at 
		         \code{dep.end}.
		smooth: does smooth on the dependent probes with smoothing factor given at \code{smooth.lambda}, finds the value of smooth 
		        for each independent probe and does the integration on them. Only needed when \code{method = "smooth"}, default \code{smooth.lambda = 2}}
		 
	\item{dep.end}{\code{numeric} or \code{character} either the name of the column \dQuote{end} in the dependent data or, when not available, an numeric value which 
	              indicates the end deviating from the start. When a numeric value is inserted, the function will do: 
	              \eqn{start + dep.end = end}. Only needed when \code{method = "window"} or \dQuote{overlap}.}
	              
	\item{window}{numeric values. Window to search for overlapping independent features per dependent probe. 
	              First value is the number of positions to the left from the middle of the probe, the second value is the number 
	              of positions to the right from the middle of the probe. Only needed when \code{method = "window"}.}
	              
	\item{smooth.lambda}{\code{numeric} factor used for smoothing the dependent data. Only needed when \code{method} = "smooth". See \link[quantsmooth]{quantsmooth} for more information.
	                     By default the \code{segment = min(nrow(dep.data), 100)}.}
	
	\item{adjust}{\code{formula} a formula like ~gender, where gender is a vector of the same size as samples. The regression models is correct for the gender effect,
	              see \link[globaltest]{gt}.}
	 
	\item{run.name}{\code{character} name of the analysis. The results will be stored in a folder with this name in the current working directory 
	                (use \code{getwd()} to print the current working directory). If missing the default folder \code{"analysis\_results"} will be generated.}
	                
	\item{...}{additional arguments for \link[globaltest]{gt} e.g. \code{model="logistic"} or when 
	           \code{permutations > 0} the null distribution is estimated using permutations, see \link[globaltest:gt]{gt}. See Details.} 
}

\details{
	The Integrated Analysis is a regression of the independent data on the dependent features. The regression itself is done using the 
	\link[globaltest]{gt}, which means that the genes in a region (e.g. a chromosome arm) are tested as a gene set. 
The individual associations between each dependent and each independent feature are calculated as Z-scores (standardized influences, see \code{\link[globaltest:gt]{?gt}}).
  
	This function splits the datasets into separate sets for each region (as specified by the \code{input.regions}) and runs the analysis for each region
separately.
  
	When running the Integrated Analysis for a predefined input region, like \dQuote{all arms} 
and \dQuote{all chrs}, output can be obtained for all input regions, as well as 
subsets of it. But note that the genomic unit must be the same: if \code{integrated.analysis} 
was run using chromosomes as units, any of the functions and plots must also use chromosomes 
as units, and not chromosome arms. Similarly, if \code{integrated analysis} was run using 
chromosome arms as units, these units must also be used to produce plots and outputs. 
For example if the \code{input.regions = "all arms"} was used, P-value plots 
(see \link{sim.plot.pvals.on.region} can be produced by inserting the \code{input.regions = "all arms"}, 
but also for instance \dQuote{1p} or \dQuote{20q}. However, to produce a plot of the whole 
chromosome, for example chromosome 1, the integrated should be re-run with \code{input.region=1}.
The same goes for \dQuote{all chrs}: P-value plots etc. can be produced for chromosome 1,2 and so on... 
but to produce plots for an arm, the \code{integrated.analysis} should be re-run for that region. 
This also goes for subregions of the chromosome like "chr1:1-1000000". 

By default the \link[globaltest]{gt} uses a \dQuote{linear} model, only when the dependent data is a \code{logical matrix} containing 
\code{TRUE} and \code{FALSE} a \dQuote{logistic} model is selected. All other models need \code{model = ""}, see \link[globaltest]{gt} 
for available models.  

}
\value{
 	No values are returned. Instead, the results of the analysis are stored
in the subdirectories of the directory specified in \code{run.name}. E.g. the z-score matrices 
are saved in subfolder \code{method}. 

	The following functions can be used to visualize the data: 
  
  \item{1) }{\link{sim.plot.zscore.heatmap} (only possible when \code{zscores = TRUE})}
  \item{2) }{\link{sim.plot.pvals.on.region}}
  \item{3) }{\link{sim.plot.pvals.on.genome}}
  \item{4) }{\link{sim.plot.overlapping.indep.dep.features}}

  \item{}{Other functions can be used to tabulate the results:}
  \item{1) }{\link{tabulate.pvals}}
  \item{2) }{\link{tabulate.top.dep.features}}
  \item{3) }{\link{tabulate.top.indep.features} (only possible when \code{zscores = TRUE}}
  \item{4) }{\link{getoverlappingregions} (only possible when tablulate.top.dep.features and tabulate.top.indep.features were run.}

}

\references{
Menezes RX, Boetzer M, Sieswerda M, van Ommen GJ, Boer JM (2009). 
Integrated analysis of DNA copy number and gene expression microarray data using gene sets.
\emph{BMC Bioinformatics}, \bold{10}, 203-.

Goeman JJ, van de Geer SA, de Kort F, van Houwelingen HC (2004). 
A global test for groups of genes: testing association with a clinical outcome. 
\emph{Bioinformatics}, \bold{20}, 93-109.
}

\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\link{SIM}, 
\link{sim.plot.zscore.heatmap}, 
\link{sim.plot.pvals.on.region}, 
\link{sim.plot.pvals.on.genome}, 
\link{tabulate.pvals}, 
\link{tabulate.top.dep.features}, 
\link{tabulate.top.indep.features}, 
\link{getoverlappingregions}, 
\link{sim.plot.overlapping.indep.dep.features},
\link[globaltest]{gt}
}

\examples{
#first run example(assemble.data)
data(samples) 
#perform integrated analysis without Z-scores using the method = "full"
integrated.analysis(samples=samples, 
					input.regions="8q", 
					zscores=FALSE, 
					method="full", 
					run.name="chr8q")
}
\keyword{multivariate}

