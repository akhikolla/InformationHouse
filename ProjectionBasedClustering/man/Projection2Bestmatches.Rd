\name{Projection2Bestmatches}
\alias{Projection2Bestmatches}

\title{
Projection to Bestmatches
}
\description{
Transformation of projected points to bestmatches defined by generalized Umatrix
}
\usage{
Projection2Bestmatches(ProjectedPoints)
}

\arguments{
  \item{ProjectedPoints}{
[1:n,1:2] n projected points in two-dimensional space.
}
}
\details{
It is assumed that an unambiguous assignment between projected points and data points is given. 

}
\value{

\item{Bestmatches}{
[1:n,1:2]   Positions of GridConverted Projected Points, which can be used for the generalized Umatrix, to the predefined Grid by Lines and Columns. First Columns has the content of the Line No and second Column of the Column number.
}
\item{LC }{[1:2] vector if Line No. and ColumnNo. which defines the size of the grid of the generalized Umatrix}

}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 
}
\author{
Michael Thrun
}
\note{
Details of the equations used are written down in [Thrun, 2018, p. 47].
}


\seealso{
\code{\link[GeneralizedUmatrix]{XYcoords2LinesColumns}}
}
\examples{
data('Hepta')
ProjList=MDS(Hepta$Data)
trafo=Projection2Bestmatches(ProjList$ProjectedPoints)
}
