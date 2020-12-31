\name{Isomap}
\alias{Isomap}
\title{Isomap}
\usage{
Isomap(Distances,k,OutputDimension=2,PlotIt=FALSE,Cls)
}
\arguments{
\item{Distances}{Symmetric [1:n,1:n] distance matrix, e.g. \code{as.matrix(dist(Data,method))}}
\item{k}{number of k nearest neighbors, if the data is fragmented choose an higher k}
 
\item{OutputDimension}{Number of dimensions in the output space, default = 2}
\item{PlotIt}{
	Default: FALSE, If TRUE: Plots the projection as a 2d visualization. 
	 If OutputDimension > 2 only the first two dimensions will be shown.
}
\item{Cls}{Optional and only relevant if PlotIt=TRUE. Numeric vector, given Classification in numbers: every element is the cluster number of a certain corresponding element of data.
}
}
\value{ProjectedPoints[1:n,OutputDimension]   n by OutputDimension matrix containing coordinates of the Projection: A matrix of the fitted configuration..}
\description{
 Isomap procetion as introduced in 2000 by Tenenbaum, de Silva and Langford

Even with a manifold structure, the sampling must be even and dense so that dissimilarities along a manifold are shorter than across the folds. If data do not have such a manifold structure, the results are very sensitive to parameter values.

}
\note{
A wrapper enabling a planar projection of the manifold learning method based on the isomap of the package vegan

if Data fragmented choose an higher k

 You can use the standard \code{ShepardScatterPlot} or the better approach through the \code{ShepardDensityPlot} of the CRAN package \code{DataVisualizations}.
}
 \author{Michael Thrun}
 \details{
 An short overview of different types of projection methods can be found in [Thrun, 2018, p.42, Fig. 4.1] (\url{https://doi.org/10.1007/978-3-658-20540-9}).
 }
\keyword{Isomap}

\concept{Manifold Learning}
\examples{
data('Hepta')
Data=Hepta$Data

Proj=Isomap(as.matrix(dist(Data)),k=7)

\dontrun{
PlotProjectedPoints(Proj$ProjectedPoints,Hepta$Cls)
}
}
\keyword{DR}