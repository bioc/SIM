\name{acgh.heatmap.with.heatmap}
\alias{acgh.heatmap.with.heatmap}
\title{Function to plot the copy number data as heatmap along side the association heatmap}
\description{
   Internal function for the SIM package that plots the association heatmap with a left panel 
   representing the log-ratios of the copy number data per sample in a color scale 
   (red is loss, green is gain). This function is used by \code{heatmap.with.plot}, which 
   in turn is used by \code{sim.plot.zscore.heatmap}.
}
\usage{
acgh.heatmap.with.heatmap(subtype, input.region, scale.heatmap, run.name)
}

\arguments{
  \item{subtype}{ Describe \code{subtype} here }
  \item{input.region}{ Describe \code{input.region} here }
  \item{scale.heatmap}{ Describe \code{scale.heatmap} here }
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
\references{ put references to the literature/web site here  }

\author{Marten Boetzer, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\note{ further notes 
 Make other sections like Warning with \\section\{Warning\}\{....\} 
}
\seealso{\code{\link{sim.plot.zscore.heatmap}}, \code{\link{heatmap.with.plot}} }
\examples{
#
}
\keyword{internal}

