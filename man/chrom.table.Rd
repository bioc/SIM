\name{chrom.table}

\alias{chrom.table}

\docType{data}

\title{Table with chromosome information}

\description{
	A table indicating the base positions of the beginning and end of chromosome arms and bands. 
   	Currently based on the UCSC March 2006/NCBI 36 build of the human genome.
}

\usage{
	data(chrom.table)
}

\format{
	A data frame with 862 observations on the following 6 variables.
  	\describe{chr, arm, band, start, end, stain}
}

\details{
	Possibly the \code{chrom.table} can be update by \link{sim.update.chrom.table}.
   	Currently only human genome support without mitochondrial DNA. 
}

\seealso{
	\link{sim.update.chrom.table}
}

\examples{
	data(chrom.table)
}

\keyword{datasets}
