\name{UniformManifoldApproximationProjection}
\alias{UniformManifoldApproximationProjection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Uniform Manifold Approximation and Projection
}
\description{
Uniform manifold approximation and projection is a technique for dimension reduction. The algorithm was described by [McInnes et al., 2018].
}
\usage{
UniformManifoldApproximationProjection(DataOrDistances, k,

Epochs,OutputDimension=2,Algorithm='umap_pkg',PlotIt=FALSE,Cls,\dots)
}

\arguments{
\item{DataOrDistances}{
Numerical matrix defined as either 

\code{Data}, i.e., [1:n,1:d], nonsymmetric, and consists of n cases of d-dimensional data points with every case having d attributes, variables or features,

or

\code{Distances}, i.e.,[1:n,1:n], symmetric and consists of n cases, e.g., \code{as.matrix(dist(Data,method))}
}
  \item{k}{
number of k nearest neighbors, Important parameter, if not given, settings of package \pkg{umap} will be used, default of package \pkg{umap} is currently 15
}
  \item{Epochs}{
Number of eppochs (scalar), i.e, training length, default of package \pkg{umap} is currently 200
}
  \item{OutputDimension}{
Number of dimensions in the Outputspace, default=2
}
	\item{Algorithm}{
	'umap': provides an interface for two implementations. One is written from scratch other one requires python \pkg{umap}
	
    'uwot_pkg': complete re-implementation in R (and C++, via the 'Rcpp' package) of \pkg{uwot}
}

  \item{PlotIt}{
	Default: FALSE, If TRUE: Plots the projection as a 2d visualization. 
	 OutputDimension>2: only the first two dimensions will be shown
}
  \item{Cls}{
Optional,: only relevant if PlotIt=TRUE. Numeric vector, given Classification in numbers: every element is the cluster number of a certain corresponding element of data.
}
  \item{\dots}{
one of the other 21 parameters that can be specified, please see \code{\link[umap]{umap.defaults}} of package \pkg{umap} for details or parameters to be set in package \pkg{uwot} depending on the choice of \code{Algorithm}.
}
}
\details{
To the knowledge of the author of this function no peer-reviewed publication of the method exists. Use with greate care.
}
\value{
List of 
 \item{ProjectedPoints}{[1:n,OutputDimension], n by OutputDimension matrix containing coordinates of the Projection}
  \item{ModelObject}{output of  \code{\link[umap]{umap}}}
   \item{Setting}{specific settings used in \code{UniformManifoldApproximationProjection}}
}
\references{
[McInnes et al., 2018]  McInnes, L., Healy, J., & Melville, J.: Umap: Uniform manifold approximation and projection for dimension reduction, arXiv preprint arXiv:1802.03426, 2018.

[Ultsch/Siemon, 1990]  Ultsch, A., & Siemon, H. P.: Kohonen's Self Organizing Feature Maps for Exploratory Data Analysis, International Neural Network Conference, pp. 305-308, Kluwer Academic Press, Paris, France, 1990.
}
\author{
Michael Thrun
}
\note{
Uniform Manifold Approximation and Projection and U-matrix [Ultsch/Siemon, 1990] are both sometimes abbreviated with Umap. Hence the abbreveviation is omitted here.
}

\seealso{
\code{\link[umap]{umap}}
}
\examples{
data('Hepta')
Data=Hepta$Data

Proj=UniformManifoldApproximationProjection(Data)

\dontrun{
PlotProjectedPoints(Proj$ProjectedPoints,Hepta$Cls)
}
}


\keyword{UniformManifoldApproximationProjection}

\concept{Uniform Manifold Approximation Projection}
\concept{Projection Method}
\concept{Dimensionality Reduction}