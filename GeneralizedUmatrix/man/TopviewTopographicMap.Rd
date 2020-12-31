\name{TopviewTopographicMap}
\alias{TopviewTopographicMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Topview of Topographic Map ind 2D
}
\description{
Fast Visualization of the Generalized U-matrix in 2D which visualizes high-dimensional distance and density based structurs of the combination two-dimensional scatter plots (projections) with high-dimensional data.
}
\usage{
TopviewTopographicMap(GeneralizedUmatrix, BestMatchingUnits,

Cls, ClsColors = NULL, Imx = NULL, Names = NULL, BmSize = 6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{GeneralizedUmatrix}{(1:Lines,1:Columns),  [1:Lines,1:Columns] Umatrix to be plotted, numerical matrix storing the U-heights, see [Thrun, 2018] for definition.}
\item{BestMatchingUnits}{(1:n,1:2), Positions of bestmatches to be plotted onto the Umatrix}
  \item{Cls}{(1:n), numerical vector of classification of \code{k} classes for the bestmatch at the given point}
  \item{ClsColors}{Vector of colors that will be used to colorize the different classes}
 \item{Imx}{a mask (Imx) that will be used to cut out the umatrix}
\item{Names}{If set: [1:k] character vector naming the k classes for the legend.
. In this case, further parameters with the possibility to adjust are:
\code{NamesCex}: (size); \code{NamesPosition}: Legend position; \code{NamesTitle}: title of legend; \code{NamesColors}: colors if ClsColors are not default (NULL).
}

\item{BmSize}{size(diameter) of the points in the visualizations. The points represent the BestMatchingUnits}

\item{\dots}{
\describe{
 \item{Tiled}{Should the Umatrix be drawn 4times?}
  \item{main}{set specific title in plot}
   \item{ExtendBorders}{scalar, extends Umatrix by toroidal continuation of the given Umatrix}
  \item{_}{Further Arguments relevant for interactive shiny application}
}

}
}
\details{
Please see \code{\link{plotTopographicMap}}. This function is currently still experimental because not all functionallity is fully tested yet.
}
\value{
plotly handler
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

[Thrun et al., 2016]  Thrun, M. C., Lerch, F., Loetsch, J., & Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision (WSCG), Vol. 24, Plzen, http://wscg.zcu.cz/wscg2016/short/A43-full.pdf, 2016.

}
\author{
Tim Schreier, Luis Winckelmann, Michael Thrun
}
\note{
Names are currently under development, Imx in testing phase.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotTopographicMap}}
}
\examples{

data("Chainlink")
Data=Chainlink$Data
Cls=Chainlink$Cls
InputDistances=as.matrix(dist(Data))
res=cmdscale(d=InputDistances, k = 2, eig = TRUE, add = FALSE, x.ret = FALSE)
ProjectedPoints=as.matrix(res$points)
#see also ProjectionBasedClustering package for other common projection methods
\donttest{
resUmatrix=GeneralizedUmatrix(Data,ProjectedPoints)
## visualization
TopviewTopographicMap(GeneralizedUmatrix = resUmatrix$Umatrix,resUmatrix$Bestmatches)
}
\dontshow{
	resUmatrix=GeneralizedUmatrix(Data[1:2,],ProjectedPoints[1:2,])
	TopviewTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches)
}
}
\keyword{topographic}% use one of  RShowDoc("KEYWORDS")
\keyword{map}% __ONLY ONE__ keyword per line
\keyword{visualization}% __ONLY ONE__ keyword per line
\keyword{U-Matrix}% __ONLY ONE__ keyword per line
\keyword{Umatrix}% __ONLY ONE__ keyword per line
\keyword{GeneralizedUmatrix}% __ONLY ONE__ keyword per line
\concept{topographic map}
\concept{Generalized U-matrix}