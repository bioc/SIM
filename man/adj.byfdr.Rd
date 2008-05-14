\name{adj.byfdr}
\alias{adj.byfdr}
\title{Adjust p-values using the Benjamini-Yekutieli method}
\description{
   Internal function for the SIM package to adjust p-values according to 
   Benjamini and Yekutieli's FDR.
}
\usage{
adj.byfdr(pvals, chrn = 0)
}

\arguments{
  \item{pvals}{ Describe \code{pvals} here }
  \item{chrn}{ Describe \code{chrn} here }
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
\references{Benjamini Y, Yekutieli D (2001). The control of the false discovery rate
in multiple testing under dependency. ANN STAT 29 (4): 1165-1188.}
\author{Marten Boetzer, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}
\note{ further notes 

 Make other sections like Warning with \section{Warning }{....} 
}
\seealso{\code{\link{sim.plot.pvals.on.region}},\code{\link{sim.plot.pvals.on.genome}}, 
\code{\link{sim.plot.adj.pvals}}, \code{\link{sim.plot.zscore.heatmap}}, 
\code{\link{tabulate.pvals}}, \code{\link{tabulate.top.dep.features}}, 
\code{\link{tabulate.top.indep.features}},\code{\link[multtest]{multtest}}  
}
\examples{
#
}
\keyword{internal}

