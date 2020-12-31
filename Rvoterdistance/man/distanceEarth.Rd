\name{distanceEarth}
\alias{distanceEarth}
\title{Calculate the distance between two points}
\description{Calculates the distance between two points on Earth, in Haversines}
\usage{
distanceEarth(lat1d, lon1d, lat2d, lon2d)
}
\arguments{
	\item{lat1d}{Numeric/Double, Latitude coordinate of point 1}
	\item{lon1d}{Numeric/Double, Longitude coordinate of point 1}
	\item{lat2d}{Numeric/Double, Latitude coordinate of point 2}
	\item{lon2d}{Numeric/Double, Longitude coordinate of point 2}
}

\value{
Haversine distance output, in meters
}
\references{
Haversine: CC Robusto, 1957
}
%\author{
%Loren Collingwood <loren.collingwood@ucr.edu>
%}
\examples{
data(king_dbox)
# Calculate distance between two points
distanceEarth(king_geo$Residence_Addresses_Latitude[1], 
king_geo$Residence_Addresses_Longitude[1], 
dbox$lat[1], dbox$long[1])

}