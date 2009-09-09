\name{sim.update.chrom.table}

\alias{sim.update.chrom.table}

\title{Update the chromomosome table}

\description{
   A function to update the genomic positions of chromosome arms. Base locations 
   of the start and end of chromosome arms should be used from the same 
   organism and build of genome as the location provided as annotation with 
   the datasets.
}
\usage{
	sim.update.chrom.table(db = "homo_sapiens_core_40_36b")
}

\arguments{
  \item{db}{database name}
}

\details{
   This functions requires \link[RMySQL]{RMySQL}. Currently \link{SIM} only supports integrated analysis on the 
   human genome without mitochondrial DNA.  
}

\value{
	Chromosome table \link{chrom.table}.
}

\references{\url{http://www.ensembl.org/info/data/mysql.html}}

\author{Marten Boetzer, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\link{SIM},
\link{chrom.table}
}

\examples{
#youn need internet connection for this!
#sim.update.chrom.table(db = "homo_sapiens_core_40_36b")
}

\keyword{database}

