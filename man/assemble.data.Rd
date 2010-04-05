\name{assemble.data}

\alias{assemble.data}

\title{Assemble the data to run the integrated analysis}

\description{
	Assembles the dependent and independent data and annotation of the both data sets.
}

\usage{
assemble.data(dep.data, 
		indep.data,		
		dep.id = "ID", 
		dep.chr = "CHROMOSOME", 
		dep.pos = "STARTPOS",
		dep.ann = NULL,		
		dep.symb,  
		indep.id = "ID", 
		indep.chr = "CHROMOSOME", 
		indep.pos = "STARTPOS",
		indep.ann = NULL,
		indep.symb,
		overwrite = FALSE, 
		run.name = "analysis_results")
}

\arguments{
	\item{dep.data}{The dependent data (\code{data.frame}), along with annotations. 
   					Each row should correspond to one feature. The following columns are expected to exist, 
   					and the column names should be inserted in the function.    
    				\code{dep.id}: A unique identifier.
    				\code{dep.chr}: The number of the chromosome (1,2, ..., 22, X, Y).
    				\code{dep.pos}: The base pair position, relative to the chromosome. 
    				\code{dep.symb}: Gene symbol (optional).
    				\code{dep.ann}: Annotation can be multiple columns.} 
    				
  \item{indep.data}{\code{data.frame} The independent data, along with annotations. 
  					Each row should correspond to one feature. The following columns are expected to exist, 
  					and the column names should be inserted in the function.    
    				\code{indep.id}: A unique identifier.
    				\code{indep.chr}: The number of the chromosome (1,2, ..., 22, X, Y).
    				\code{indep.pos}: The base pair position, relative to the chromosome. 
    				\code{indep.symb}: Gene symbol (optional).
    				\code{indep.ann}: Annotation can be multiple columns.}
    				    
  \item{dep.ann}{\code{vector} with either the names of the columns or the 
                 column numbers in the dependent data that contain the annotation.}

  \item{indep.ann}{\code{vector} with either the names of the columns or the 
                   column numbers in the independent data that contain the annotation.}

  \item{dep.id}{\code{vector} with the column name in the dependent 
                data that contains the ID. Will be used in the \link{sim.plot.zscore.heatmap} function. 
                Empty ID's will be substituted by \code{NA}.}

  \item{dep.chr}{\code{vector} with column name in the dependent 
                 data that contains the chromosome numbers.}

  \item{dep.pos}{\code{vector} with the column name in the dependent data 
                 that contains the position on the chromosome in bases.}

  \item{dep.symb}{Optional, either missing or a single vector with the column name 
                  in the dependent data that contains the symbols. Will be used in \link{sim.plot.zscore.heatmap} 
                  as label.}

  \item{indep.id}{\code{vector} with the column name in the independent 
                  data that contains the ID. Will be used in the \link{sim.plot.zscore.heatmap} function. 
                  Empty ID's will be substituted by \code{NA}.}

 \item{indep.chr}{\code{vector} with the column name in the independent data 
                  that contains the chromosome numbers.}

  \item{indep.pos}{\code{vector} with the column name in the independent data 
                   that contains the position on the chromosome in bases.}

  \item{indep.symb}{Optional, either missing or a vector with the column name 
                    in the dependent data that contains the Symbols. Will be used in \link{sim.plot.zscore.heatmap} 
                    as label.}

  \item{overwrite}{\code{logical}, indicate when a \code{run.name} is already present, 
                    the results can be overwritten.}

  \item{run.name}{Name of the analysis. The results will be 
                  stored in a folder with this name in the current working directory 
                  (use \code{getwd()} to print the current working directory). 
                   If the missing, the default folder \code{"analysis_results"} will be generated.}
}

\details{
	Based on the chromosome and probe position an absolute position is calculated according to 
	\eqn{chromosome number * 1e9 + probe position}. Chromosome column is converted to \code{factor} and releveled according to
	the levels of the \link{chrom.table}, so the only levels allowed are \code{c(1:22, "X", "Y")}. 
	Currently only human genome support without mitochondrial DNA. 
}

\value{
	No values are returned. Instead, the datasets and annotation columns are stored in 
separate files in the \code{data} folder in the directory specified in \code{run.name}. 
If \code{assemble.data} has run successfully, the \link{integrated.analysis} function can be performed.
}
\author{Marten Boetzer, Melle Sieswerda, Renee X. de Menezes  \email{R.X.Menezes@lumc.nl}}

\seealso{
\link{SIM}, \link{integrated.analysis} 
}

\examples{
#load the datasets and the samples to run the integrated analysis
data(expr.data)
data(acgh.data)
data(samples) 
         
#assemble the data
assemble.data(dep.data = acgh.data, 
              indep.data = expr.data,              
              dep.id="ID", 
              dep.chr = "CHROMOSOME",
              dep.pos = "STARTPOS",
              dep.ann = colnames(acgh.data)[1:4],    
              dep.symb="Symbol",                          
              indep.id="ID",
              indep.chr = "CHROMOSOME", 
              indep.pos = "STARTPOS",
              indep.ann = colnames(expr.data)[1:4], 
              indep.symb="Symbol",                
              overwrite = TRUE,
              run.name = "chr8q")
}
\keyword{manip}

