\name{pmatrixForEsom}
\alias{pmatrixForEsom}
\title{pmatrixForEsom}
\arguments{
\item{Data}{A \code{[n,k]} matrix containing the data}
\item{Weights}{Weights stored as a list in a 2D matrix}
\item{Lines}{Number of lines of the SOM that is described by weights}
\item{Columns}{Number of columns of the SOM that is described by weights}
\item{Radius}{The radius for measuring the density within the hypersphere}
\item{PlotIt}{If set the Pmatrix will also be plotted}
\item{Toroid}{Are BestMatches placed on a toroid grid? TRUE by default}
}
\value{UstarMatrix}
\description{Generates a Pmatrix based on the weights of an ESOM.}

\references{
Ultsch, A.: Maps for the visualization of high-dimensional data spaces, Proc. Workshop on Self organizing Maps (WSOM), pp. 225-230, Kyushu, Japan, 2003.

Ultsch, A., Loetsch, J.: Computed ABC Analysis for Rational Selection of Most Informative Variables in Multivariate Data, PloS one, Vol. 10(6), pp. e0129767. doi 10.1371/journal.pone.0129767, 2015.

Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016.

}

\examples{
data("Hepta")
e = esomTrain(Hepta$Data, Key = 1:nrow(Hepta$Data))
Pmatrix = pmatrixForEsom(Hepta$Data,
                         e$Weights,
                         e$Lines,
                         e$Columns,
                         e$Toroid)
plotMatrix(Pmatrix, ColorStyle = "Pmatrix")
}

\keyword{Pmatrix}
\keyword{P-matrix}
\keyword{P-Matrix}
