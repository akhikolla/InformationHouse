\name{ListAsEsomNeurons}
\alias{ListAsEsomNeurons}
\title{Converts List to WTS}
\arguments{
\item{wts_list[1:(Lines*Columns),1:Variables]}{Matrix with weights in the 2nd dimension(not list() like in R)}
\item{Lines}{Lines/Height of the desired grid}
\item{Columns}{Columns/Width of the desired grid}
}
\value{
\item{EsomNeurons[1:Lines, 1:Columns, 1:Variables]}{3 dimensional array containing the weights of the neural grid. For a more general explanation see reference}}
\description{Converts wts data in list form into a 3 dimensional array}
\details{One could describe this function as a transformation or a special case of long to wide format, see also \code{\link{EsomNeuronsAsList}}}
\author{Michael Thrun, Florian Lerch}
\references{
Ultsch, A.: Maps for the visualization of high-dimensional data spaces, Proc. Workshop on Self organizing Maps (WSOM), pp. 225-230, Kyushu, Japan, 2003. }