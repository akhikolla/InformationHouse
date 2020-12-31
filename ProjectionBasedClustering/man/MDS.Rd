\name{MDS}
\alias{MDS}
\title{Multidimensional Scaling (MDS)}
\usage{
MDS(DataOrDistances,method='euclidean',OutputDimension=2,PlotIt=FALSE,Cls)
}

\arguments{
\item{DataOrDistances}{
Numerical matrix defined as either 

\code{Data}, i.e., [1:n,1:d], nonsymmetric, and consists of n cases of d-dimensional data points with every case having d attributes, variables or features,

or

\code{Distances}, i.e.,[1:n,1:n], symmetric and consists of n cases, e.g., \code{as.matrix(dist(Data,method))}
}
\item{method}{method specified by distance string: 'euclidean','cityblock=manhatten','cosine','chebychev','jaccard','minkowski','manhattan','binary'}
\item{OutputDimension}{Number of dimensions in the Outputspace, default=2}
  \item{PlotIt}{
	Default: FALSE, If TRUE: Plots the projection as a 2d visualization. 
	}
	\item{Cls}{[1:n,1] Optional,: only relevant if PlotIt=TRUE. Numeric vector, given Classification in numbers: every element is the cluster number of a certain corresponding element of data.
}
}
\value{

  \item{ProjectedPoints}{
 [1:n,OutputDimension], n by OutputDimension matrix containing coordinates of the Projectio
}

  \item{Eigenvalues}{
 the eigenvalues of MDSvalues*MDSvalues'
}
  \item{Stress}{
 Shephard-Kruskal Stress  
}
}
\description{
 Classical multidimensional scaling of a data matrix. Also known as principal coordinates analysis}
\author{Michael Thrun}
\note{
A wrapper for \code{cmdscale}

 You can use the standard \code{ShepardScatterPlot} or the better approach through the \code{ShepardDensityPlot} of the CRAN package \code{DataVisualizations}.
}
 \details{
 An short overview of different types of projection methods can be found in [Thrun, 2018, p.42, Fig. 4.1] (\url{https://doi.org/10.1007/978-3-658-20540-9}).
 }
 
  \examples{
data('Hepta')
Data=Hepta$Data

Proj=MDS(Data)

\dontrun{
PlotProjectedPoints(Proj$ProjectedPoints,Hepta$Cls)
}
}

\keyword{MDS}

\concept{Classical multidimensional scaling}
\concept{Projection Method}
\concept{Dimensionality Reduction}
\keyword{DR}