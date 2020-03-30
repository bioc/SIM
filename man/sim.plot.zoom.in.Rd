\name{sim.plot.zoom.in}

\alias{sim.plot.zoom.in}

\title{Zoom in on heatmap}

\description{
	Zoom in on pervious produced heatmap by \link{sim.plot.zscore.heatmap}
}

\usage{
sim.plot.zoom.in(call) 
}

\arguments{
 \item{call}{\code{language}, function call of \link{sim.plot.zscore.heatmap}}
}

\details{
	\link{sim.plot.zscore.heatmap} returns (invisible) the function call and as attribute the local function environment. 
	So by adjusting the the xlim and ylim a zoom in is created.
}

\value{an additional plot is created, on the previous plot a rectangle indicating the zoomed region.}

\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes \email{R.X.Menezes@lumc.nl}}

\seealso{
\link{SIM}, 
\link{sim.plot.zscore.heatmap}
}

\examples{
#first run example(assemble.data)
#and example(integrated.analysis)

#plot the zscores in a heatmap
sim.plot <- sim.plot.zscore.heatmap(input.regions = "8q",  
                                    adjust.method = "BY", 
                                    run.name = "chr8q", pdf = FALSE)
#only when runned interactive
if(interactive())
	sim.plot.zoom.in(sim.plot)

}
\keyword{iplot}


