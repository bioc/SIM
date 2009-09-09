\name{acgh.data}

\alias{acgh.data}

\docType{data}

\title{Array Comparative Genomic Hybridization data}

\description{
	Copy number data taken from Pollack et al. PNAS. 2002, 99(20): 12963-8.
}

\usage{
	data(acgh.data)
}

\format{
	A data frame with 99 observations on 45 variables. The first 4 columns are the unique identifier, symbol
  	for the chromosome and start position of the probe the next 41 columns are the copy number measurements of 41 samples.  
}

\details{
	A subset of the original data is taken namely all the probes on the long arm of chromosome 8.	 
}

\source{
	\url{http://genome-www.stanford.edu/aCGH_breast/data.shtml} 
}

\references{
Pollack JR, Sorlie T, Perou CM, Rees CA, Jeffrey SS, Lonning PE, Tibshirani R, Botstein D, Borresen-Dale AL, Brown PO (2002).
Microarray analysis reveals a major direct role of DNA copy number alteration in the transcriptional program of human breast tumors.
\emph{Proc Natl Acad Sci USA}. \bold{1}, 99(20), 12963-8. 	 
}

\examples{
	data(acgh.data)
}

\keyword{datasets}