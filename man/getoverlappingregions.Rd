\name{getoverlappingregions}

\alias{getoverlappingregions}

\title{Get the overlapping regions between independent and dependent regions}

\description{
	Generates a table with overlapping regions.
}

\usage{
getoverlappingregions(independent_regions, 
                      dependent_regions, 
                      method = c("union", "overlapping"))
}

\arguments{
	\item{independent_regions}{data.frame().Independent regions found with tab \link{tabulate.top.indep.features}.}

	\item{dependent_regions}{data.frame().Independent regions found with tab \link{tabulate.top.dep.features}.}

	\item{method}{method to estimate the overlapping regions, either "union" or "overlapping".
                  \code{overlapping} outputs only the overlapping parts of the overlapping regions.
                  \code{union} outputs the whole region. Say independent region = 1-10, dependent region = 5-12. The output is 1-12.}
}

\details{
Calculates the overlap between two tables.
}

\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}


\seealso{
\link{SIM}, 
\link{tabulate.top.dep.features}, 
\link{tabulate.top.indep.features}, 
\link{sim.plot.overlapping.indep.dep.features}
}

\examples{
#no examples yet!
}

\keyword{misc}

