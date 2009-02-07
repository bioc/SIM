\name{quantsmooth.heatmap}
\alias{quantsmooth.heatmap}

\title{Function to plot smoothed copy number data along side association heatmap}
\description{
   Internal function for the SIM package used by \code{heatmap_with_plot}, which 
   in turn is used by \code{sim.plot.zscore.heatmap}.
}
\usage{
quantsmooth.heatmap(input.region, dep.data, lambda, run.name)
}

\arguments{
  \item{input.region}{ Describe \code{input.region} here }
  \item{dep.data}{ Describe \code{dep.data} here }
  \item{lambda}{ Describe \code{lambda} here }
  \item{run.name}{ Describe \code{run.name} here }
}
\details{
   See the help in the main function sim.plot.zscore.heatmap. 
}
\value{
  Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}

}
\references{Jan Oosting, Paul Eilers and Renee Menezes (2006). quantsmooth: Quantile smoothing and
  genomic visualization of array data. R package version 1.2.0.}
\author{Marten Boetzer, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}
\note{ further notes 

 Make other sections like Warning with \section{Warning }{....}  
}
\seealso{\code{\link{heatmap.with.plot}},\code{\link{sim.plot.zscore.heatmap}},\code{\link[quantsmooth]{quantsmooth}} }
\examples{
#
}
\keyword{internal}

