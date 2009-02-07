\name{get.abs.start}
\alias{get.abs.start}

\title{get.abs.start}
\description{
   Internal function for the SIM package usesd by the \code{read.data} function. 
}
\usage{
get.abs.start(chr, pos, dataset)
}

\arguments{
  \item{chr}{ Describe \code{chr} here }
  \item{pos}{ Describe \code{pos} here }
  \item{dataset}{ Describe \code{dataset} here }
}
\details{
   This function generates unique chromosome positions by taking the start position 
   of a probe in base pairs, multiplying by 1,000,000,000 and adding the chromosome number.
}
\value{
  Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}

}
\references{ put references to the literature/web site here  }
\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\note{ further notes 

 Make other sections like Warning with \\section\{Warning\}\{...\} 
}
\seealso{ objects to See Also as \code{\link{help}}}
\examples{
#
}
\keyword{internal}

