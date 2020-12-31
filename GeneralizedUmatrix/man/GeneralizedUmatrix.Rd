\name{GeneralizedUmatrix}
\alias{GeneralizedUmatrix}
\alias{getUmatrix4Projection}
\title{Generalized U-Matrix for Projection Methods published in [Thrun/Ultsch, 2020]}
\usage{
GeneralizedUmatrix(Data,ProjectedPoints,

PlotIt=FALSE,Cls=NULL,Toroid=TRUE,Tiled=FALSE,ComputeInR=FALSE)
}
\arguments{
\item{Data}{[1:n,1:d] array of data: n cases in rows, d variables in columns}
\item{ProjectedPoints}{[1:n,2] matrix containing coordinates of the Projection: A matrix of the fitted configuration.}
\item{PlotIt}{
 Optional,bool, defaut=FALSE, if =TRUE: U-Marix of every current Position of Databots will be shown
 However, the amount of details shown will be less than in \code{\link{plotTopographicMap}}.
}
\item{Cls}{
Optional, For plotting, see \code{plotUmatrix} in package Umatrix
}
\item{Toroid}{
Optional, Default=TRUE, 

==FALSE planar computation with borders defined by projection method

==TRUE: toroid borderless (toroidal) computation, the four borders defined by projection method are ignored.
}
\item{Tiled}{
Optional,For plotting see \code{plotUmatrix} in package Umatrix
}
\item{ComputeInR}{
Optional,  =T: Rcode, =F Cpp Code
}

}
\value{List with
\item{Umatrix}{
 [1:Lines,1:Columns]  [1:Lines,1:Columns] Umatrix to be plotted, numerical matrix storing the U-heights, see [Thrun, 2018] for definition.
}
\item{EsomNeurons}{
[1:Lines,1:Columns,1:weights]  3-dimensional numeric array (wide format), not wts (long format)
}
\item{Bestmatches}{
[1:n,1:2]   Positions of GridConverted Projected Points on the Umatrix to the predefined Grid by Lines and Columns, First Columns has the content of the Line No and second Column of the Column number.
}
\item{gplotres}{
 Ausgabe von ggplot
}

}
\description{Generalized U-Matrix visualizes high-dimensional distance and density based structurs in two-dimensional scatter plots of projectios methods like CCA, MDS, PCA or NeRV [Ultsch/Thrun, 2017] with the help of a topographic map with hypsometrioc tints [Thrun et al. 2016] using a simplified emergent SOM published in [Thrun/Ultsch, 2020]. 

}

\author{Michael Thrun}
\details{
Introduced first in the PhD thesis in [Thrun, 2018, p.46]. Furthermore the two parts of the work were peer-reviewed and published in [Ultsch/Thrun, 2017, Thrun/Ultsch, 2020].
}
\references{
[Thrun et al., 2016]  Thrun, M. C., Lerch, F., Loetsch, J., & Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision (WSCG), Vol. 24, Plzen, http://wscg.zcu.cz/wscg2016/short/A43-full.pdf, 2016.

[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

[Ultsch/Thrun, 2017]  Ultsch, A., & Thrun, M. C.: Credible Visualizations for Planar Projections, in Cottrell, M. (Ed.), 12th International Workshop on Self-Organizing Maps and Learning Vector Quantization, Clustering and Data Visualization (WSOM), IEEE Xplore, France, 2017.

[Thrun/Ultsch, 2020]  Thrun, M. C., & Ultsch, A.: Uncovering High-Dimensional Structures of Projections from Dimensionality Reduction Methods, MethodsX, Vol. in press, pp. 101093. doi 10.1016/j.mex.2020.101093, 2020.

}
\examples{
data("Chainlink")
Data=Chainlink$Data
Cls=Chainlink$Cls
InputDistances=as.matrix(dist(Data))
res=cmdscale(d=InputDistances, k = 2, eig = TRUE, add = FALSE, x.ret = FALSE)
ProjectedPoints=as.matrix(res$points)
\dontrun{
Stress = ProjectionBasedClustering::KruskalStress(InputDistances,

as.matrix(dist(ProjectedPoints)))
}
\dontshow{
	resUmatrix=GeneralizedUmatrix(Data[1:2,],ProjectedPoints[1:2,])
	#plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches)#testing takes longer than 5 secs
}
\donttest{
	resUmatrix=GeneralizedUmatrix(Data,ProjectedPoints)
	plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,Cls)
}
}
\keyword{visualization}% __ONLY ONE__ keyword per line
\keyword{U-Matrix}% __ONLY ONE__ keyword per line
\keyword{Umatrix}% __ONLY ONE__ keyword per line
\keyword{GeneralizedUmatrix}% __ONLY ONE__ keyword per line
\concept{Generalized U-matrix}