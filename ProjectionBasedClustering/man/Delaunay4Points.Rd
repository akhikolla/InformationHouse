\name{Delaunay4Points}
\alias{Delaunay4Points}
\title{
Adjacency matrix of the delaunay graph for BestMatches of Points
}
\description{
Calculates the adjacency matrix of the delaunay graph for BestMatches (BMs) in tiled form if BestMatches are located on a toroid grid
}
\usage{
Delaunay4Points(Points, IsToroid = TRUE,Grid=NULL,PlotIt=FALSE)
}

\arguments{
  \item{Points}{[1:n,1:3] matrix containing the BMKey, X and Y coordinates of the n,  BestMatches NEED NOT BE UNIQUE,  however, there is an edge in the Deaunay between duplicate points!}
   \item{IsToroid}{OPTIONAL, logical, indicating if BM's are on a toroid grid. Default is True}
    \item{Grid}{OPTIONAL, A vector of length 2, containing the number of lines and columns of the Grid}
     \item{PlotIt}{OPTIONAL, bool, Plots the graph}
}
\details{
as
}
\value{
 Delaunay[1:n,1:n]               adjacency matrix of the Delaunay-Graph
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, ISBN: 978-3-658-20539-3, Heidelberg, 2018.
}
\author{
Michael Thrun
}


\keyword{Delaunay}
\keyword{Points}% __ONLY ONE__ keyword per line
