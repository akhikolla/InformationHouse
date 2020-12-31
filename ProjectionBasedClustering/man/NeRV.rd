\name{NeRV}
\alias{NeRV}
\title{
Neighbor Retrieval Visualizer (NeRV) 
}
\description{
Projection is done by the neighbor
retrieval visualizer (NeRV) 
}
\usage{
NeRV(Data, lambda = 0.1, neighbors = 20, iterations = 10, 

cg_steps = 2, cg_steps_final = 40, randominit = T, OutputDimension = 2,

PlotIt = FALSE, Cls)
}
\arguments{
  \item{Data}{Numerical matrix of the Data to be projected, [1:n,1:d], nonsymmetric, and consists of n cases of d-dimensional data points with every case having d attributes, variables or features}
  \item{lambda}{Optional: Controls the trustworthiness-continuity tradeoff. Default = 0.1}
  \item{neighbors}{Optional: Set the number of nearest neighbours that each point should have. Should be positive. Default = 20}
  \item{iterations}{Optional: The number of iterations to perform. Default = 10}
  \item{cg_steps}{Optional: The number of conjugate gradient steps to perform per iteration in NeRV's optimization scheme. Default = 2}
  \item{cg_steps_final}{Optional: The number of conjugate gradient steps to perform on the final iteration in NeRV's optimization scheme. Default = 40}
  \item{randominit}{Optional: TRUE: Random Initialization (default), FALSE: PCA initializiation}
  \item{OutputDimension}{Optional: Number of dimensions in the Outputspace, default=2}
  \item{PlotIt}{Optional: Should the projected points be plotted? Default: FALSE. Note: this is only usefull if OutputDimension = 2.}
  \item{Cls}{Optional: Vector containing the number of the class for each row in Data. This is only used to color the points according to their classes if PlotIt = T}
	}


\details{
Uses the NeRV projection with matrix Data and lambda. Lambda controls the trustworthiness-continuity tradeoff.

 An short overview of different types of projection methods can be found in [Thrun, 2018, p.42, Fig. 4.1] (\url{https://doi.org/10.1007/978-3-658-20540-9}).

}
\value{
OutputDimension-dimensional matrix of projected points
}
\note{
PCA initialization changes form the original C++ Sourcecode of \url{http://research.cs.aalto.fi/pml/software/dredviz/} to the R version of the projections package.
Other changes are made only regarding data types of Rcpp in comparison to the original C++ Source code.

 You can use the standard \code{ShepardScatterPlot} or the better approach through the \code{ShepardDensityPlot} of the CRAN package \code{DataVisualizations}.
}
\references{
Jarkko Venna, Jaakko Peltonen, Kristian Nybo, Helena Aidos, and Samuel Kaski. Information Retrieval Perspective to Nonlinear Dimensionality Reduction for Data Visualization. Journal of Machine Learning Research, 11:451-490, 2010.

Jarkko Venna and Samuel Kaski. Nonlinear Dimensionality Reduction as Information Retrieval. In Marina Meila and Xiaotong Shen, editors, Proceedings of AISTATS 2007, the 11th International Conference on Artificial Intelligence and Statistics. Omnipress, 2007. JMLR Workshop and Conference Proceedings, Volume 2: AISTATS 2007.
}
\author{
Michael Thrun, Felix Pape
}

\examples{
data('Hepta')
Data=Hepta$Data
\dontrun{
Proj=NeRV(Data)
PlotProjectedPoints(Proj$ProjectedPoints,Hepta$Cls)
}

\dontshow{
NeRV(Data[1:10,])
}
}

\keyword{NeRV}

\concept{neighbor retrieval visualizer}
\concept{Projection Method}
\concept{Dimensionality Reduction}
\keyword{DR}