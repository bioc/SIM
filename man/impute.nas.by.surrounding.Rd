\name{impute.nas.by.surrounding}

\alias{impute.nas.by.surrounding}

\title{Impute NA's in array-CGH data}

\description{
	Replace an NA by the median of the surrounding features in the same sample.
}

\usage{
impute.nas.by.surrounding(dataset, 
                          window.size = 5)
}

\arguments{
	\item{dataset}{\link{data.frame} with the dataset to replace the NA's by the medians of the surrounding
                 features.}
                 
	\item{window.size}{numeric value, specifying of how many features around the NA the median should be taken.}
}

\details{
	This function can be used when the dependent dataset in the \link{integrated.analysis} 
function is array-CGH data and contains probes that have an NA. To avoid loosing data by 
throwing away the probes with NA's, the impute.nas.by.surrounding function can be 
used which simply takes the median of the probes around an NA. The number of probes used 
for the imputation is chosen by giving a value for window.size. This script takes quite long to run!}

\value{
A \code{data.frame} is returned, containing the inserted \dQuote{dataset} all \code{NA} replaced with median of the window of 
size \dQuote{window.size} around the \code{NA}.
}

\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\link{SIM}, 
\link{assemble.data}, 
\link{integrated.analysis} 
}

\examples{
#no examples yet!
}

\keyword{manip}
