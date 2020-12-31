\name{umatrixForEsom}
\alias{umatrixForEsom}
\title{umatrixForEsom}
\arguments{
	\item{Weights}{Weights from which the Umatrix will be calculated}
	\item{Lines}{Number of lines of the SOM that is described by weights}
	\item{Columns}{Number of columns of the SOM that is described by weights}
	\item{Toroid}{Boolean describing if the neural grid should be borderless}
}
\value{Umatrix}
\description{Calculate the Umatrix for given ESOM projection}

\examples{
data("Hepta")
e = esomTrain(Hepta$Data, Key = 1:nrow(Hepta$Data))
umatrix = umatrixForEsom(e$Weights,
                         Lines=e$Lines,
                         Columns=e$Columns,
                         Toroid=e$Toroid)
plotMatrix(umatrix,e$BestMatches)
}

\references{Ultsch, A. and H.P. Siemon, Kohonen's Self Organizing Feature Maps for Exploratory Data Analysis. 1990.}
\keyword{Umatrix}
\keyword{U-matrix}
\keyword{U-Matrix}
