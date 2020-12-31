\name{PolarSwarm}
\alias{PolarSwarm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polar Swarm (Pswarm)
}
\description{
Swarm-based Projection method using game theory published in [Thrun/Ultsch, 2020].
}
\usage{
PolarSwarm(DataOrDistances, method = "euclidean", PlotIt = FALSE, Cls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{DataOrDistances}{
Numerical matrix defined as either 

\code{Data}, i.e., [1:n,1:d], nonsymmetric, and consists of n cases of d-dimensional data points with every case having d attributes, variables or features,

or

\code{Distances}, i.e.,[1:n,1:n], symmetric and consists of n cases, e.g., \code{as.matrix(parallelDist::parallelDist(Data,method))}
}
  \item{method}{
	If \code{Data} is given the method to computing the distances can be specified here.
	Please see the documentation of package \pkg{parallelDist} for the types that are possible.
}
  \item{PlotIt}{
	Default: FALSE, If TRUE: Plots the projection as a 2d visualization. 
	 OutputDimension>2: only the first two dimensions will be shown
}
  \item{Cls}{
Optional,: only relevant if PlotIt=TRUE. Numeric vector, given Classification in numbers: every element is the cluster number of a certain corresponding element of data.
}
}
\details{
By exploting swarm intelligence and game theory no parameter have to be set.
}
\value{
List of 
 \item{ProjectedPoints}{[1:n,2], n by 2 matrix containing coordinates of the Projection}
  \item{ModelObject}{output of  \code{\link[DatabionicSwarm]{Pswarm}}}
}
\references{
[Thrun/Ultsch, 2020]  Thrun, M. C., & Ultsch, A.: Swarm Intelligence for Self-Organized Clustering, Artificial intelligence, Vol. in press, pp. doi 10.1016/j.artint.2020.103237, 2020.
}
\author{
Michael Thrun
}


\seealso{
\code{\link[DatabionicSwarm]{Pswarm}}
}
\examples{
data('Hepta')
Data=Hepta$Data
\donttest{
Distances=as.matrix(dist(Data))
Proj=PolarSwarm(Data)}
\dontrun{
PlotProjectedPoints(Proj$ProjectedPoints,Hepta$Cls)}
\dontshow{Proj=PolarSwarm(as.matrix(dist(Data[1:20,])))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{Pswarm}
\keyword{swarm}
\keyword{DR}
\keyword{PolarSwarm}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{Polar Swarm}
\concept{Projection Method}
\concept{Dimensionality Reduction}
% Use only one concept per line.
