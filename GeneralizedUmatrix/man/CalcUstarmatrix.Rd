\name{CalcUstarmatrix}
\alias{CalcUstarmatrix}
\alias{calcUstarmatrix}
\title{ Calculate the U*matrix for a given Umatrix and Pmatrix. }
\arguments{
\item{Umatrix[1:Lines,1:Column])}{ Local averages of distances at each point of the trainedGridWts[1:Lines,1:Column,1:variables] of ESOM or other SOM of same format}
\item{Pmatrix[1:Lines,1:Column]}{ Local densities at each point of the trainedGridWts[1:Lines,1:Column,1:variables] of ESOM or other SOM of same format}
}
\value{UStarMatrix[1:Lines,1:Column]}
\description{
 Calculate the U*matrix for a given Umatrix and Pmatrix.}
\author{Michael Thrun}
\references{Ultsch, A. U* C: Self-organized Clustering with Emergent Feature Maps. in Lernen, Wissensentdeckung und Adaptivitaet (LWA). 2005. Saarbruecken, Germany.}
\keyword{U*matrix}
\keyword{U*-matrix}
\keyword{U*-Matrix}
\keyword{Ustarmatrix}
\keyword{Ustar-matrix}
\keyword{Ustar-Matrix}