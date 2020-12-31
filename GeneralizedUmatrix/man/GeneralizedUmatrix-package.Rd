\name{GeneralizedUmatrix-package}
\alias{GeneralizedUmatrix-package}
\docType{package}
\title{
\packageTitle{GeneralizedUmatrix}
}
\description{
\packageDescription{GeneralizedUmatrix}
}
\details{
For a brief introduction to \pkg{GeneralizedUmatrix} please see the vignette \href{../doc/GeneralizedUmatrix.html}{Introduction of the Generalized Umatrix Package}.

For further details regarding the generalized Umatrix see [Thrun, 2018], chapter 4-5, or [Thrun/Ultsch, 2020].

If you want to verifiy your clustering result externally, you can use \code{Heatmap} or \code{SilhouettePlot} of the CRAN package \code{DataVisualizations}.

\packageIndices{GeneralizedUmatrix}
}
\author{
Michal Thrun

Maintainer: Michael Thrun <mthrun@informatik.uni-marburg.de>
}
\references{
[Thrun/Ultsch, 2020]  Thrun, M. C., & Ultsch, A.: Uncovering High-Dimensional Structures of Projections from Dimensionality Reduction Methods, MethodsX, Vol. in press, pp. 101093. doi 10.1016/j.mex.2020.101093, 2020.

[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

[Ultsch/Thrun, 2017]  Ultsch, A., & Thrun, M. C.: Credible Visualizations for Planar Projections, in Cottrell, M. (Ed.), 12th International Workshop on Self-Organizing Maps and Learning Vector Quantization, Clustering and Data Visualization (WSOM), IEEE Xplore, France, 2017.
}
\examples{
data("Chainlink")
Data=Chainlink$Data
Cls=Chainlink$Cls
InputDistances=as.matrix(dist(Data))
res=cmdscale(d=InputDistances, k = 2, eig = TRUE, add = FALSE, x.ret = FALSE)
ProjectedPoints=as.matrix(res$points)
#see also ProjectionBasedClustering package for other common projection methods
#see DatabionicSwarm for projection method without parameters or objective function
# ProjectedPoints=DatabionicSwarm::Pswarm(Data)$ProjectedPoints
\donttest{
resUmatrix=GeneralizedUmatrix(Data,ProjectedPoints)
plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,Cls)
}
##Interactive Island Generation 
## from a tiled Umatrix (toroidal assumption)
\dontrun{
	Imx = ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(resUmatrix$Umatrix,
	resUmatrix$Bestmatches)
	plotTopographicMap(resUmatrix$Umatrix,

	resUmatrix$Bestmatches, Imx = Imx)
}
#External Verification
\dontrun{

 DataVisualizations::Heatmap(Data,Cls)
 #if spherical cluster strcuture
 DataVisualizations::SilhouettePlot(Data,Cls)
}

}
\keyword{GeneralizedUmatrix}
\keyword{U-matrix}
\keyword{Umatrix}
\keyword{U-Matrix}
\keyword{Visualization}
\keyword{DimensionalityReduction}
\keyword{DR}
\keyword{Projection}