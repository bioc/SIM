\name{link.metadata}

\alias{link.metadata}

\title{Link a metadata annotation file to expression ID}

\description{
	Get annotation out of a AnnotationData package and link them to the expression data 
using the expression probe ID's
}

\usage{
link.metadata(data = expr.data, 
              col.ID.link = 1, 
              chr = as.list(hgu133plus2CHR), 
              chrloc = as.list(hgu133plus2CHRLOC), 
              symbol = as.list(hgu133plus2SYMBOL))
}

\arguments{ 
	\item{data}{\code{data.frame} with expression data including an expression probe ID column.}
	
    \item{col.ID.link}{numeric value, specifying the column of \code{data} that contains the ID to link with the \code{poslist}.}
    
    \item{chr}{\code{list} specifying the metadata annotation of the chromosome location on the genome.}
    
    \item{chrloc}{\code{list} specifying the metadata annotation of the  location of the probe on the chromosome.}
    
    \item{symbol}{\code{list} specifying the metadata annotation of the symbol corresponding to the probe.}
}
\details{
Often, the annotation for expression array probes lack chromosome position information. Therefore, 
this function adds the two required columns to run the \link{integrated.analysis}:
"CHROMOSOME" and "STARTPOS". In addition, the optional column, "Symbol" is added.
}

\value{A \code{data.frame} is returned, containing a dataset with annotation columns 
which can be used for\link{integrated.analysis}. 
}

\author{Marten Boetzer, Melle Sieswerda, Renee x Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\link{RESOURCERER.annotation.to.ID}
}

\examples{

# first download and install the AnnotationData package for your expression array platform
# for example
\dontrun{library(hgu133plus2)}
\dontrun{expr.data <- link.metadata(data, col.ID.link = 1, chr = as.list(hgu133plus2CHR), 
chrloc = as.list(hgu133plus2CHRLOC), symbol = as.list(hgu133plus2SYMBOL))}
}
\keyword{database}

