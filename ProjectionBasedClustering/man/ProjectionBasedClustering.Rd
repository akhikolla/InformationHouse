\name{ProjectionBasedClustering}
\alias{ProjectionBasedClustering}
\title{
Automatic Projection-based Clustering (PBC) [Thrun/Ultsch, 2020]
}
\description{
Three steps are necessary for PBC. First, a projection method has to be chosen to generate projected points of high-dimensional data points. Second, the generalized U*-matrix has to be applied to the projected points by using a simplified emergent self-organizing map (ESOM) algorithm which is an unsupervised neural network [Thrun, 2018]. The resulting generalized U-matrix can be visualized by the topographic map [Thrun et al., 2016]. Third, the clustering itself is built on top of the generalized U-matrix using the concept of the abstract U-Matrix and shortest graph paths using \code{ShortestGraphPathsC}. 
}
\usage{
ProjectionBasedClustering(k, DataOrDistances, BestMatches, LC,

StructureType = TRUE, PlotIt = FALSE, method = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{number of clusters, how many to you see in the 3d landscape?}
  \item{DataOrDistances}{Numerical matrix that will be used for clustering with one DataPoint per row, defined as either as

\code{Data}, i.e., [1:n,1:d], nonsymmetric, and consists of n cases of d-dimensional data points with every case having d attributes, variables or features,

or

\code{Distances}, i.e.,[1:n,1:n], symmetric and consists of n cases, e.g., \code{as.matrix(dist(Data,method))}

}
  \item{BestMatches}{[1:n,1:2] Matrix with positions of Bestmatches=ProjectedPoints, one matrix line per data point}
  \item{LC}{grid size c(Lines,Columns)}
  \item{StructureType}{Optional, bool; =TRUE: compact structure of clusters assumed, =FALSE: connected structure of clusters assumed. For the two options vor Clusters, see [Thrun, 2017] or Handl et al. 2006}
  \item{PlotIt}{Optional, bool, Plots Dendrogramm}
  \item{method}{Optional, distance method used in \pkg{parallelDist} if \code{Data} given.}
}
\details{
ProjectionBasedClustering is a flexible and robust clustering framework based on a chose projection method and
 projection method a parameter-free high-dimensional data visualization technique. The visualization combines projected points with a topographic map with hypsometric colors, defined by the generalized U-matrix  (see package GeneralizedUmatrix function GeneralizedUmatrix). 
 
 The clustering method with no sensitive parameters is done in this function and the algorithm is introduced in detail in [Thrun/Ultsch, 2020]. The clustering can be verified by the visualization and vice versa. If you want to verifiy your clustering result externally, you can use \code{Heatmap} or \code{SilhouettePlot} of the CRAN package \code{DataVisualizations}.
 
If \pkg{parallelDist} is not installed, function automatically falls back to \code{\link{dist}}.
}
\value{
Cls                 [1:n] vector with selected classes of the bestmatches. You can use \code{plotTopographicMap(Umatrix,Bestmatches,Cls)} for verification.
}
\examples{
data('Hepta')
#Step I: 2d projection
\donttest{projectionpoints=NeRV(Hepta$Data)

#Step II (Optional): Computation of Generalized Umatrix
library(GeneralizedUmatrix)
visualization=GeneralizedUmatrix(Data = Hepta$Data,projectionpoints)
# Visualizuation of GeneralizedUmatrix
library(GeneralizedUmatrix)
TopviewTopographicMap(visualization$Umatrix,visualization$Bestmatches)
#or in 3D if rgl package exists
#plotTopographicMap(visualization$Umatrix,visualization$Bestmatches)

# Step III: Automatic Clustering
trafo=Projection2Bestmatches(projectionpoints)
# number of Cluster from dendrogram  (PlotIt=T) or visualization above
Cls=ProjectionBasedClustering(k=7, Hepta$Data, 

trafo$Bestmatches, trafo$LC,PlotIt=TRUE)

# Verification of Clustering
TopviewTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls)
#or in 3D if rgl package exists
#plotTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls)
}
}
\references{
[Thrun et al., 2016]  Thrun, M. C., Lerch, F., Lötsch, J., & Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision (WSCG), Vol. 24, Plzen, http://wscg.zcu.cz/wscg2016/short/A43-full.pdf, 2016.

[Thrun/Ultsch, 2017] Thrun, M.C., Ultsch, A.: Projection based Clustering, Conf. Int. Federation of Classification Societies (IFCS),DOI:10.13140/RG.2.2.13124.53124, Tokyo, 2017.

[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

[Thrun/Ultsch, 2020]  Thrun, M. C., & Ultsch, A.: Using Projection based Clustering to Find Distance and Density based Clusters in High-Dimensional Data, Journal of Classification, Vol. in press, Springer, DOI: 10.1007/s00357-020-09373-2, 2020. 
}
\author{
Michael Thrun
}
\note{
Often it is better to mark the outliers  manually after the prozess of clustering; use in this case the visualization \code{plotTopographicMap} of the package GeneralizedUmatrix. If you would like to mark the outliers interactivly in the visualization use 
the \code{interactiveClustering} function.
}
\keyword{PBC}
\concept{Projection-based Clustering}
\keyword{clustering}% __ONLY ONE__ keyword per line
\concept{cluster analysis}% __ONLY ONE__ keyword per line
\keyword{cluster}% __ONLY ONE__ keyword per line