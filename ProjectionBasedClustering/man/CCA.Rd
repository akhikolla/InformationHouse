\name{CCA}
\alias{CCA}
\title{Curvilinear Component Analysis (CCA)}
\usage{
CCA(DataOrDistances,Epochs,OutputDimension=2,method='euclidean',

alpha0 = 0.5, lambda0,PlotIt=FALSE,Cls)
}
\arguments{
\item{DataOrDistances}{
Numerical matrix defined as either 

\code{Data}, i.e., [1:n,1:d], nonsymmetric, and consists of n cases of d-dimensional data points with every case having d attributes, variables or features,

or

\code{Distances}, i.e.,[1:n,1:n], symmetric and consists of n cases, e.g., \code{as.matrix(dist(Data,method))}
}
	\item{Epochs}{Number of eppochs (scalar), i.e, training length}
	\item{OutputDimension}{Number of dimensions in the Outputspace, default=2}
	\item{method}{method specified by distance string. One of: 'euclidean','cityblock=manhatten','cosine','chebychev','jaccard','minkowski','manhattan','binary'}
	\item{alpha0}{(scalar) initial step size, 0.5 by default}
	\item{lambda0}{(scalar) initial radius of influence, 3*max(std(D)) by default}

	\item{PlotIt}{
	Default: FALSE, If TRUE: Plots the projection as a 2d visualization. 
	 OutputDimension>2: only the first two dimensions will be shown
	}
	\item{Cls}{[1:n,1] Optional,: only relevant if PlotIt=TRUE. Numeric vector, given Classification in numbers: every element is the cluster number of a certain corresponding element of data.
}
}
\value{A n by OutputDimension matrix containing coordinates of the projected points.}
\description{
  CCA Projects data vectors using Curvilinear Component Analysis [Demartines/Herault, 1995],[Demartines/Herault, 1997].
 
   Unknown values (NaN's) in the data: projections of vectors with
 unknown components tend to drift towards the center of the
  projection distribution. Projections of totally unknown vectors are
  set to unknown (NaN).
	}
 \author{Florian Lerch}
 
 \note{
 Only Transfered from matlab to R. Matlabversion: Contributed to SOM Toolbox 2.0, February 2nd, 2000 by Juha Vesanto
  \url{http://www.cis.hut.fi/projects/somtoolbox/}

 You can use the standard \code{ShepardScatterPlot} or the better approach through the \code{ShepardDensityPlot} of the CRAN package \code{DataVisualizations}.
}
 \details{
 An short overview of different types of projection methods can be found in [Thrun, 2018, p.42, Fig. 4.1] (\url{https://doi.org/10.1007/978-3-658-20540-9}).
 }
 \references{
[Demartines/Herault, 1997]  Demartines, P., & Herault, J.: Curvilinear component analysis: A self-organizing neural network for nonlinear mapping of data sets, IEEE Transactions on Neural Networks, Vol. 8(1), pp. 148-154. 1997.
	
[Demartines/Herault, 1995]  Demartines, P., & Herault, J.: CCA:" Curvilinear component analysis", Proc. 15 Colloque sur le traitement du signal et des images, Vol. 199, GRETSI, Groupe d'Etudes du Traitement du Signal et des Images, France 18-21 September, 1995.


 }
\examples{
data('Hepta')
Data=Hepta$Data

Proj=CCA(Data,Epochs=20)

\dontrun{
PlotProjectedPoints(Proj$ProjectedPoints,Hepta$Cls)
}
}
\keyword{CCA}

\concept{Curvilinear Component Analysis}
\concept{Projection Method}
\concept{Dimensionality Reduction}
\keyword{DR}