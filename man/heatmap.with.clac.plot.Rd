\name{heatmap.with.clac.plot}
\alias{heatmap.with.clac.plot}
\title{heatmap.with.clac.plot}
\description{
   Internal function for the SIM package used by the sim.plot.zscore.heatmap 
   when run with plot.method = \code{"clac"}. 
}
\usage{
heatmap.with.clac.plot(input.region, dep.data, Normal.data, windowsize = windowsize, run.name)
}

\arguments{
  \item{input.region}{ Describe \code{input.region} here }
  \item{dep.data}{ Describe \code{dep.data} here }
  \item{Normal.data}{ Describe \code{Normal.data} here }
  \item{windowsize}{ Describe \code{windowsize} here }
  \item{run.name}{ Describe \code{run.name} here }
}
\details{
   We used code from the clac package to implement this function.
}
\value{
  Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}

}
\references{ Pei Wang and with contributions from Balasubramanian Narasimhan.
  (2004). clac: Clust Along Chromosomes, a method to call gains/losses
  in CGH array data. R package version 0.1-2.
  }
\author{Marten Boetzer, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\note{ further notes 

 Make other sections like Warning with \\section\{Warning\}\{...\} 
}
\seealso{ objects to See Also as \code{\link{help}}}
\examples{
#
}
\keyword{internal}

