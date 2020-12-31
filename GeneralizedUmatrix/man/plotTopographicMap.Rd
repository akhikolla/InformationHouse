\name{plotTopographicMap}
\alias{plotTopographicMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualizes the Generalized U-matrix in 3D
}
\description{
The generalized U-matrix is visualized as the topographic map with hypsometric tints. The topographic map represents high-dimensional distance and density-based structurs in form of a 3D landscape.
}
\usage{
plotTopographicMap(GeneralizedUmatrix, BestMatchingUnits,

Cls=NULL,ClsColors=NULL,Imx=NULL,Names=NULL,

BmSize=0.5,RenderingContourLines=TRUE,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{GeneralizedUmatrix}{(1:Lines,1:Columns),  [1:Lines,1:Columns] Umatrix to be plotted, numerical matrix storing the U-heights, see [Thrun, 2018] for definition.}
\item{BestMatchingUnits}{(1:n,1:2), Positions of bestmatches to be plotted as spheres onto the topographic map}
  \item{Cls}{(1:n), numerical vector of classification of \code{k} clusters, one label for each bestmatch at that given point}
  \item{ClsColors}{Vector of colors that will be used to colorize the different clusters, default is GeneralizedUmatrix::DefaultColorSequence}
 \item{Imx}{a mask (Imx) that will be used to cut out the umatrix}
\item{Names}{If set: [1:k] character vector naming the k clusters for the legend.
. In this case, further parameters with the possibility to adjust are:
\code{NamesCex}: (size); \code{NamesPosition}: Legend position; \code{NamesTitle}: title of legend; \code{NamesColors}: colors if ClsColors are not default (NULL), etc.
}

\item{BmSize}{size(diameter) of the points in the visualizations. The points represent the BestMatchingUnits}
 \item{RenderingContourLines}{FALSE: disables plotting of contour lines resulting in a much faster plot.}
\item{\dots}{

 Besides the legend/names parameter the list of further parameters, use only of you know what you are doing:
\describe{
 \item{Tiled}{Should the Umatrix be drawn 4times?}
 \item{ShowAxis}{shall the axis be shown?}
 \item{NoLevels}{number of contour lines}
  \item{ExtendBorders}{scalar, extends Umatrix by toroidal continuation of the given Umatrix}
 \item{Colormap}{in the case of density p matrix...}
\item{title}{same as main}
   \item{main}{same as title}
      \item{sub}{same as in \code{\link{plot}}}
         \item{xlab}{same as in \code{\link{plot}}}   
         \item{ylab}{same as in \code{\link{plot}}}
         \item{zlab}{same as in \code{\link{plot}}}
   \item{NamesPosition}{same as in \code{\link[rgl]{bgplot3d}}}
    \item{NamesColors}{same as col in \code{\link[rgl]{bgplot3d}}}
     \item{NamesCex}{same as cex in \code{\link[rgl]{bgplot3d}}}
      \item{NamesTitle}{same as title in \code{\link[rgl]{bgplot3d}}}
      \item{NamesPch}{same as pch in \code{\link[rgl]{bgplot3d}}}
 }
 
}
}
\details{
The visualization of this function is a topographic map with hypsometric tints (Thrun, Lerch, L?tsch, & Ultsch, 2016).
"Hypsometric tints are surface colors that represent ranges of elevation (Patterson and Kelso 2004). Here, contour lines are combined with a specific color scale. The color scale is chosen to display various valleys, ridges, and basins: blue colors indicate small distances (sea level), green and brown colors indicate middle distances (low hills), and white colors indicate vast distances (high mountains covered with snow and ice). Valleys and basins represent clusters, and the watersheds of hills and mountains represent the borders between clusters. In this 3D landscape, the borders of the visualization are cyclically connected with a periodicity (L,C). The number of clusters can be estimated by the number of valleys of the visualization. The clustering is valid if mountains do not partition clusters indicated by colored points of the same color and colored regions of points (see examples in section 4.1 and 4.2)."[Thrun/Ultsch, 2020].  

A central problem in clustering is the correct estimation of the number of clusters. This is addressed by the topographic map which allows assessing the number of clusters as the number of valleys (Thrun et al., 2016). Please see chapter 5 of [Thrun, 2018] for further details.
}

\value{
An object of class "htmlwidget" in mode invisible, please \code{\link[rgl]{rglwidget}} for details.
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

[Thrun et al., 2016]  Thrun, M. C., Lerch, F., Loetsch, J., & Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision (WSCG), Vol. 24, Plzen, http://wscg.zcu.cz/wscg2016/short/A43-full.pdf, 2016.

[Thrun/Ultsch, 2020]  Thrun, M. C., & Ultsch, A. : Using Projection based Clustering to Find Distance and Density based Clusters in High-Dimensional Data, Journal of Classification, DOI 10.1007/s00357-020-09373-2, in press, Springer, 2020. 
}
\author{
Michael Thrun
}
\note{
First version of algorithm was partly based on the Umatrix package.
}

\seealso{
\code{\link{GeneralizedUmatrix}}}

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
plotTopographicMap(GeneralizedUmatrix = resUmatrix$Umatrix,resUmatrix$Bestmatches)
}

## Open window in specific resolution
#relevant if Names given
\donttest{
library(rgl)
r3dDefaults$windowRect = c(0,0,1200,1200) 
plotTopographicMap(GeneralizedUmatrix = resUmatrix$Umatrix,resUmatrix$Bestmatches)
}
\dontrun{
## To save as STL for 3D printing
 rgl::writeSTL("GenerelizedUmatrix_3d_model.stl")

## Save the visualization as a picture with
library(rgl)
rgl.snapshot('test.png')
}

## Save interactive html file
\dontrun{
widgets=plotTopographicMap(GeneralizedUmatrix = resUmatrix$Umatrix,resUmatrix$Bestmatches)
if(requireNamespace("htmlwidgets"))
  htmlwidgets::saveWidget(widgets,file = "interactiveTopographicMap.html")
}
\dontshow{
	resUmatrix=GeneralizedUmatrix(Data[1:2,],ProjectedPoints[1:2,])
	#plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches)#testing takes longer than 5 secs
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{topographic}% use one of  RShowDoc("KEYWORDS")
\keyword{map}% __ONLY ONE__ keyword per line
\keyword{visualization}% __ONLY ONE__ keyword per line
\keyword{U-Matrix}% __ONLY ONE__ keyword per line
\keyword{Umatrix}% __ONLY ONE__ keyword per line
\keyword{GeneralizedUmatrix}% __ONLY ONE__ keyword per line
\concept{topographic map}
\concept{Generalized U-matrix}