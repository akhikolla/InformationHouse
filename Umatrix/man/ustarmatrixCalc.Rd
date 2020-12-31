\name{ustarmatrixCalc}
\alias{ustarmatrixCalc}
\title{ustarmatrixCalc}
\arguments{
\item{Umatrix}{A given Umatrix}
\item{Pmatrix}{A density matrix}
}
\value{UStarMatrix}
\description{The UStarMatrix is a combination of the Umatrix (average distance to neighbours) and Pmatrix (density in a point). It can be used to improve the Umatrix, if the dataset contains density based structures.}

\references{Ultsch, A. U* C: Self-organized Clustering with Emergent Feature Maps. in Lernen, Wissensentdeckung und Adaptivitaet (LWA). 2005. Saarbruecken, Germany.}

\examples{
data("Hepta")
e = esomTrain(Hepta$Data, Key = 1:nrow(Hepta$Data))
Pmatrix = pmatrixForEsom(Hepta$Data,
                         e$Weights,
                         e$Lines,
                         e$Columns,
                          e$Toroid)
Ustarmatrix = ustarmatrixCalc(e$Umatrix, Pmatrix)
plotMatrix(Ustarmatrix, e$BestMatches)
}

\keyword{U*matrix}
\keyword{U*-matrix}
\keyword{U*-Matrix}
\keyword{Ustarmatrix}
\keyword{Ustar-matrix}
\keyword{Ustar-Matrix}
