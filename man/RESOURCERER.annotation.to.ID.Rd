\name{RESOURCERER.annotation.to.ID}

\alias{RESOURCERER.annotation.to.ID}

\title{Link RESOURCERER annotation file to expression ID}

\description{Get annotation out of the RESOURCERER annotation file and link them to expression data with help of expression ID's
}

\usage{
RESOURCERER.annotation.to.ID(data, poslist, col.ID.link = 1, col.poslist.link = 1)
}

\arguments{
	\item{data}{\link{data.frame} with expression data including an expression ID column.}
	
	\item{poslist}{\link{data.frame} containing the RESOURCERER annotation file}
	
	\item{col.ID.link}{numeric value, specifying the column of \code{data} that contains the ID to link with the \code{poslist}.}
	
	\item{col.poslist.link}{numeric value, specifying the column of \code{poslist} that contains the ID to link with the \code{data}.}
}

\details{
This function will output the inserted dataset, including the necessary, for \link{integrated.analysis}, a
nnotation columns: \code{"CHROMOSOME"}, \code{"STARTPOS"}and \code{"Symbol"} out of the inserted RESOURCE
RER annotation file \code{poslist}.
}

\value{A \code{data.frame} is returned, containing a dataset with annotation columns which can be used for
 \link{integrated.analysis}}

\author{Marten Boetzer, Melle Sieswerda, Renee x Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\link{link.metadata}
}


\examples{

# download expression array annotation from RESOURCERER 
# ftp://occams.dfci.harvard.edu/pub/bio/tgi/data/Resourcerer
# it may be necessary to remove the first row, which states the genome build used for mapping
\dontrun{read.an <- read.delim("affy_U133Plus2.txt", sep="\t", header=T)} 

# get physical mapping columns
\dontrun{expr.data <- RESOURCERER_annotation_to_ID(data = read.expr, poslist = read.an, col.ID.link = 1, col.poslist.link = 1)}
}
\keyword{database}

