\name{GeneratePmatrix}
\alias{GeneratePmatrix}
\title{Generates the P-matrix}
\arguments{
\item{Data}{[1:n,1:d], A \code{[n,d]} matrix containing the data}
\item{EsomNeurons}{[1:Lines,Columns,1:Weights] 3D array of weights given by ESOM or sESOM algorithm.}
\item{Radius}{The radius for measuring the density within the hypersphere.}
\item{PlotIt}{If set the Pmatrix will also be plotted}
\item{\dots}{If set the Pmatrix will also be plotted}
}
\value{PMatrix[1:Lines,1:Columns]}
\description{Generates a P-matrix too visualize only density based structures of high-dimensional data.}
\details{
To set the Radius the ABCanalysis of high-dimensional distances can be used [Ultsch/Lötsch, 2015]. For a deteailed definition and equation of automated density estimation (Radius) see Thrun et al. 2016.
}
\author{Michael Thrun}
\references{
Ultsch, A.: Maps for the visualization of high-dimensional data spaces, Proc. Workshop on Self organizing Maps (WSOM), pp. 225-230, Kyushu, Japan, 2003. 

Ultsch, A., Loetsch, J.: Computed ABC Analysis for Rational Selection of Most Informative Variables in Multivariate Data, PloS one, Vol. 10(6), pp. e0129767. doi 10.1371/journal.pone.0129767, 2015.

Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016. 

}
\keyword{Pmatrix}
\keyword{P-matrix}
\keyword{P-Matrix}