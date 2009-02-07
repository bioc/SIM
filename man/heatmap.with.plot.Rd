\name{heatmap.with.plot}
\alias{heatmap.with.plot}

\title{heatmap.with.plot}
\description{
   Internal function for the SIM package used in the function \code{sim.plot.zscore.heatmap}. 
}
\usage{
heatmap.with.plot(dep.data, input.region = input.region, plot.method = plot.method, Normal.data = Normal.data, subtype = subtype, acgh.heatmap.scale = acgh.heatmap.scale, run.name = run.name, lambda = lambda, windowsize = windowsize, x = x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, distfun = dist, hclustfun = hclust, reorderfun = function(d, w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, verbose = getOption("verbose"), ...)
}

\arguments{
  \item{dep.data}{ Describe \code{dep.data} here }
  \item{input.region}{ Describe \code{input.region} here }
  \item{plot.method}{ Describe \code{plot.method} here }
  \item{Normal.data}{ Describe \code{Normal.data} here }
  \item{subtype}{ Describe \code{subtype} here }
  \item{acgh.heatmap.scale}{ Describe \code{acgh.heatmap.scale} here }
  \item{run.name}{ Describe \code{run.name} here }
  \item{lambda}{ Describe \code{lambda} here }
  \item{windowsize}{ Describe \code{windowsize} here }
  \item{x}{ Describe \code{x} here }
  \item{Rowv}{ Describe \code{Rowv} here }
  \item{Colv}{ Describe \code{Colv} here }
  \item{distfun}{ Describe \code{distfun} here }
  \item{hclustfun}{ Describe \code{hclustfun} here }
  \item{reorderfun}{ Describe \code{reorderfun} here }
  \item{add.expr}{ Describe \code{add.expr} here }
  \item{symm}{ Describe \code{symm} here }
  \item{revC}{ Describe \code{revC} here }
  \item{scale}{ Describe \code{scale} here }
  \item{na.rm}{ Describe \code{na.rm} here }
  \item{margins}{ Describe \code{margins} here }
  \item{ColSideColors}{ Describe \code{ColSideColors} here }
  \item{RowSideColors}{ Describe \code{RowSideColors} here }
  \item{cexRow}{ Describe \code{cexRow} here }
  \item{cexCol}{ Describe \code{cexCol} here }
  \item{labRow}{ Describe \code{labRow} here }
  \item{labCol}{ Describe \code{labCol} here }
  \item{main}{ Describe \code{main} here }
  \item{xlab}{ Describe \code{xlab} here }
  \item{ylab}{ Describe \code{ylab} here }
  \item{keep.dendro}{ Describe \code{keep.dendro} here }
  \item{verbose}{ Describe \code{verbose} here }
  \item{\dots}{ Describe \code{\dots} here }
}
\details{
   If necessary, more details than the description above 
}
\value{
  Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
}
\references{ put references to the literature/web site here  }
\author{Marten Boetzer, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}
\note{ further notes 

 Make other sections like Warning with \\section\{Warning\}\{...\} 
}
\seealso{ objects to See Also as \code{\link{sim.plot.zscore.heatmap}}}
\examples{
#
}
\keyword{internal}

