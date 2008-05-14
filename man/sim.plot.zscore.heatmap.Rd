\name{sim.plot.zscore.heatmap}
\alias{sim.plot.zscore.heatmap}
\title{Association heatmap from z-scores}
\description{Produces an association heatmap that shows the association (standardized influence) of 
 each independent feature (expression measurement) with each dependent feature 
 (copy number measurement). A p-value bar on the left indicates test signficance. 
 A color bar on top indicates genes with mean z-scores across the signficant copy 
 number probes above a set threshold. A summary of the copy number data helps to identify 
 what copy number alterations are present in a region of association with expression. 
 Positive association can mean copy number gain and increased expression, or deletion and 
 decreased expression. The heatmaps can also be used in an exploratory analysis, 
 looking for very local effects of copy number changes (usually small amplifications) on 
 gene expression, that do not lead to a significant test result.
}
\usage{
sim.plot.zscore.heatmap(input.regions = "all chrs", significance = 0.2, z.threshold = 3, show.names.indep = FALSE, show.names.dep = FALSE, adjust.method = c("BY", "BH", "raw"), scale = "auto", plot.method = c("none"), Normal.data = if (plot.method == "clac") FALSE, windowsize = 5, lambda = 2, subtype = FALSE, acgh.heatmap.scale = "auto", pdf = TRUE, run.name = NULL, ...)
}

\arguments{

  \item{input.regions}{\code{\link{vector}} indicating the regions to be analyzed. Can be defined in four ways:
  \code{1) predefined input region: } insert a predefined input region, choices are:
  \code{"all chrs"}, \code{"all chrs auto"}, \code{"all arms"}, \code{"all arms auto"}
  In the predefined regions \code{"all arms"} and \code{"all arms auto"} the arms 13p, 
  14p, 15p, 21p and 22p are left out, because in most studies there are no or few probes 
  in these regions. To include them, just make your own \code{\link{vector}} of arms.
  \code{2) whole chromosome(s): }insert a single chromosome or a list of chromosomes as a \code{\link{vector}:}
  \code{c(1, 2, 3)}.
  \code{3) chromosome arms: } insert a single chromosome arm  or a list of chromosome arms like 
  \code{c("1q", "2p", "2q")}.
  \code{4) subregions of a chromosome: } insert a chromosome number followed by the start and end position
 like \code{c("chr1_1-1000000")}
  These regions can also be combined, e.g. \code{c("chr1_1-1000000","2q", 3)}.
  See \code{details} for more information.}
  
  \item{significance}{Threshold to select the significant dependent
   features. Only these features are used to calculate the mean z-scores per independent 
   feature (expression probe).}

  \item{z.threshold}{Threshold to display a green or red bar in the color bar on top of 
  the heatmap for independent features with mean z-scores above \code{z.threshold} (high positive 
  association) or below \code{-z.threshold} (high negative association).}

  \item{show.names.indep}{\code{Boolean}. If set to TRUE, displays the names (\code{indep.id} and \code{in
dep.symb} entered in the \code{\link{assemble.data}}) of the independent features with mean z-scores above or 
below the \code{z.threshold} in the heatmap.}

  \item{show.names.dep}{\code{Boolean}. If set to TRUE, displays the names (\code{dep.id} and \code{dep.sy
mb} entered in the \code{\link{assemble.data}}) of the \code{significant} dependent features in the heatmap.}

  \item{adjust.method}{Method used to adjust the p-values for multiple
   testing. Either \code{"BY"} (recommended when copy number is used as dependent data), 
   \code{"BH"} or \code{"raw"}. See \code{\link{SIM}} for more information about adjusting 
   p-values. Defaults to "BY".}
   
  \item{scale}{Vector specifying the color scale in the heatmap.
    If scale="auto", the maximum and minimum value of all z-scores will be calculated 
    and set as the limits for all analyzed regions. Another option is to define a custom scale, 
    e.g. scale = c(-5,5).}
    
  \item{plot.method}{Summary plot of copy number data in left panel. Either \code{"clac"}, 
  \code{"smooth"},\code{"heatmap"}, or \code{"none"}. Should only be used when the \code{dep.data} 
  is array-CGH. The \code{"clac"} plot is a consensus of the aberration frequencies across 
  all samples. CLAC requires at least three normal/diploid arrays. For more details see 
  \code{?clac.preparenormal.R}. The \code{"smooth"} plot smoothes the copy number log ratios 
  per sample, see \code{?quantsmooth} for more details. The \code{"heatmap"} method produces 
  an aCGH heatmap where green indicates gain, and red loss. The scale of the aCGH heatmap 
  is automatically set to the min and max of the aCGH measurements of the analyzed regions. 
  Default is plot.method = \code{"none"}, no additional plot will be drawn.}
   
  \item{Normal.data}{\code{\link{vector}}, required for plot.method = \code{"clac"}, indicating 
  least three normal samples in the dependent data. Insert the column names of 
  the samples that are normal e.g. for the first three \code{samples}: \code{Normal.data = 1:3}.
  If no normal samples are available, use Normal.data = FALSE. Then Normal.data 
  are generated by calculating probe medians of three subsets of the dependent data.}

  \item{windowsize}{Numeric value, specifying the window size to carry out the average 
  smooth for \code{plot.method="clac"}. For more details see \code{?clac.preparenormal.R}.}

  \item{lambda}{Numeric value, specifying the quantile smoothing parameter for 
  \code{plot.method="smooth"}. See \code{?quantsmooth} and \code{references} for more information.}

  \item{subtype}{This variable must be a vector with the same length as \code{samples} or FALSE. The 
  vector will be transformed to a factor and the levels of this will be coloured according to 
  their subtype. When \code{subtype}=FALSE, all the samples will be coloured black.}

  \item{acgh.heatmap.scale}{Vector specifiing the color scale in the aCGH heatmap.
    If scale="auto", the maximum and minimum value of all aCGH values will be calculated 
    and set as the limits for all analyzed regions. Another option is to define a custom scale, 
    e.g. scale = c(-5,5).}
  \item{pdf}{Logical. Indicate whether to generate a pdf of the plots in the heatmap\_zscores subdirectory or plot to screen.} 
  \item{run.name}{Name of the analysis. The results will be 
   stored in a folder with this name in the current working directory 
   (use \code{getwd()} to print the current working directory). 
   If the \code{run.name = NULL}, the default folder \code{"analysis_results"} will be generated.} 
\item{...}{additional arguments passed on to \link{image}}
}

\details{
The \code{sim.plot.zscore.heatmap} function can only run after the \code{\link{integrated.analysis}} 
is run with \code{zscores = TRUE}.

The results are returned as a single-page pdf containing an association heatmap of the regions 
listed in \code{input.regions}. For high-density arrays large files will be produced, both 
demanding more memory available from your computer to produce them as well as being heavier to 
open on screen. To avoid this, analyze chromosome arms as units instead of chromosomes, both 
here and in \code{input.regions = "all arms"}. 
  
The heatmap contains the z-scores generated by the function \code{\link{integrated.analysis}} with 
\code{zscores=T}. The dependent features are plotted from bottom to top, the independent 
features from left to right. Positive associations are shown in green, negative associations in red 
(color scale on the right). At the left side of the heatmap a color bar represents the 
multiple testing corrected p-values of the probes in the dependent data (copy number), also  
with a color legend. Dependening on which \code{plot.method} is used, a summary of copy number 
changes is shown on the left. At the top of the heatmap is a color bar corresponding to 
the mean z-scores of the independent features (expression data) that are above or below 
the \code{z.threshold}. If \code{show.names.indep} is set to TRUE, labels will be drawn for 
the probes with mean z-scores greater than \code{z.threshold} or lower than \code{-z.threshold} 
at the bottom of the heatmap. If \code{show.names.dep} is set to TRUE, labels will be drawn for 
the significant dependent probes lower than \code{significance} to the right of the heatmap.

}
\value{No values are returned. The results are stored in a subdirectory of \code{run.name} as pdf.}

\references{\item{1}{Eilers PH, de Menezes RX. Quantile smoothing of array CGH data. Bioinformatics. 2005 Apr 1;21(7):1146-53.} 
\item{2}{Wang P, Kim Y, Pollack J, Narasimhan B, Tibshirani R. A method for calling gains and losses in array CGH data. Biostatistics. 2005; 6:45-58.}
}

\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\code{\link{SIM}}, \code{\link{assemble.data}}, \code{\link{integrated.analysis}}, 
\code{\link{sim.plot.pvals.on.region}}, \code{\link{sim.plot.pvals.on.genome}}, 
\code{\link{tabulate.pvals}}, \code{\link{tabulate.top.dep.features}}, 
\code{\link{tabulate.top.indep.features}}, \code{\link{impute.nas.by.surrounding}}, 
\code{\link{sim.update.chrom.table}},\code{\link[fields]{image.plot}},\code{\link[marray]{maPalette}}
}


\examples{
#load the datasets and the samples to run the integration for
data(expr.data)
data(acgh.data)
data(samples) 
         
#assemble the data
assemble.data(dep.data = acgh.data, indep.data = expr.data, ann.dep = colnames(acgh.data)[1:4], ann.indep = colnames(expr.data)[1:4], dep.id="ID",dep.chr = "CHROMOSOME",dep.pos = "STARTPOS",dep.symb="Symbol",  indep.id="ID",indep.chr = "CHROMOSOME", indep.pos = "STARTPOS", indep.symb="Symbol", overwrite = TRUE,run.name ="chr8")

#run the integrated analysis
integrated.analysis(samples = samples, input.regions = 8, adjust=FALSE, zscores=TRUE, method = "auto", run.name= "chr8")

# use functions to plot the results of the integrated analysis

#plot the zscores in a heatmap
sim.plot.zscore.heatmap(input.regions = 8, significance=0.2, z.threshold=3, show.names.dep=TRUE, show.names.indep=TRUE, adjust.method = c("BY"), scale="auto", plot.method = "smooth", pdf = FALSE, run.name = "chr8")
}
\keyword{misc}

