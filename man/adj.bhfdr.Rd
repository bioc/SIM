\name{adj.bhfdr}
\alias{adj.bhfdr}
\title{Adjust p-values using the Benjamini-Hochberg method}
\description{
   Internal function to adjust p-values according to Benjamini and Hochberg's FDR. 
}
\usage{
adj.bhfdr(pvals, chrn = 0)
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
\references{Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate -
a practical and powerful approach to multiple testing. J Roy Stat Soc B
57 (1): 289-300.}
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

