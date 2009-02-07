\name{make.Normal.data}
\alias{make.Normal.data}

\title{Function to generate normal data for CLAC}
\description{
   This function is used when the the copy number data should be represented in 
   a CLAC plot, while there are no normal samples in the copy number data. The 
   (tumor) samples are divided into three parts, the probe medians are taken  
   for each group, thereby generating three fake normal samples. Can give 
   distortions of real copy number when aberrations occur in more than half of the samples. 
}
\usage{
make.Normal.data(data)
}

\arguments{
  \item{data}{ Describe \code{data} here }
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
\author{Marten Boetzer, Renee X. de Menezes \email{R.X.Menezes@lumc.nl}}

\note{ further notes 

 Make other sections like Warning with \\section\{Warning\}\{...\}
}
\seealso{\code{\link{heatmap.with.clac.plot}}}
\examples{
#
}
\keyword{internal}

