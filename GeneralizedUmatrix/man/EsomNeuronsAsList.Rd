\name{EsomNeuronsAsList}
\alias{EsomNeuronsAsList}
\title{Converts wts data (EsomNeurons) into the list form}
\arguments{
\item{EsomNeurons[1:Lines, 1:Columns, 1:Variables]}{high dimensional array with grid positions in the first two dimensions}
}
\value{\item{TrainedNeurons[1:(Lines*Columns),1:Variables]}{List of Weights as a matrix (not \code{\link{list}} like in R) as matrix or two dimensional array}}
\details{One could describe this function as a transformation or a special case of wide to long format, see also \code{\link{ListAsEsomNeurons}}}
\description{Converts wts data into the list form}
\author{Michael Thrun, Florian Lerch}
\references{
Ultsch, A. Maps for the visualization of high-dimensional data spaces. in Proc. Workshop on Self organizing Maps. 2003.
}