\name{esomTrain}
\alias{esomTrain}
\title{Train an ESOM (emergent self organizing map) and project data}

\arguments{
	\item{Data}{Data that will be used for training and projection}
	\item{Lines}{Height of grid}
	\item{Columns}{Width of grid}
	\item{Epochs}{Number of Epochs the ESOM will run}
	\item{Toroid}{If TRUE, the grid will be toroid}
	\item{NeighbourhoodFunction}{Type of Neighbourhood; Possible values are: "cone", "mexicanhat" and "gauss"}
	\item{StartLearningRate}{Initial value for LearningRate}
	\item{EndLearningRate}{Final value for LearningRate}
	\item{StartRadius}{Start value for the Radius in which will be searched for neighbours}
	\item{EndRadius}{End value for the Radius in which will be searched for neighbours}
	\item{NeighbourhoodCooling}{Cooling method for radius; "linear" is the only available option at the moment}
	\item{LearningRateCooling}{Cooling method for LearningRate; "linear" is the only available option at the moment}
	\item{shinyProgress}{Generate progress output for shiny if Progress Object is given}
	\item{ShiftToHighestDensity}{If True, the Umatrix will be shifted so that the point with highest density will be at the center}
	\item{InitMethod}{name of the method that will be used to choose initializations
	Valid Inputs: "uni_min_max": uniform distribution with minimum and maximum from sampleData
	"norm_mean_2std": normal distribuation based on mean and standard deviation of sampleData}
	\item{Key}{Vector of numeric keys matching the datapoints. Will be added to Bestmatches}
	\item{UmatrixForEsom}{If TRUE, Umatrix based on resulting ESOM is calculated and returned}
}

\value{List with
	\item{BestMatches}{BestMatches of datapoints}
	\item{Weights}{Trained weights}
	\item{Lines}{Height of grid}
	\item{Columns}{Width of grid}
	\item{Toroid}{TRUE if grid is a toroid}
}

\description{The ESOM (emergent self organizing map) algorithm as defined by [Ultsch 1999]. A set of weigths(neurons) on a two-dimensional grid get trained to adapt the given datastructure. The weights will be used to project data on a two-dimensional space, by seeking the BestMatches for every datapoint.}

\details{On a toroid grid, opposing borders are connected.}

\references{
Kohonen, T., Self-organized formation of topologically correct feature maps. Biological cybernetics, 1982. 43(1): p. 59-69.

Ultsch, A., Data mining and knowledge discovery with emergent self-organizing feature maps for multivariate time series. Kohonen maps, 1999. 46: p. 33-46.
}

\examples{
data('Hepta')
res=esomTrain(Hepta$Data, Key = 1:nrow(Hepta$Data))
}

\keyword{ESOM}
\keyword{esom}
\keyword{som}
\keyword{SOM}
\keyword{self organizing map}
