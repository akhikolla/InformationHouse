
\name{king_geo}
\alias{king_geo}
\title{Dataset of King County voters' lat/long }
\description{Dataset of 5,000 randomly selected King County voters' lat/long, as of 2016}
\usage{
data(king_dbox)
}
\format{A dataset of 5,000 rows and two columns:
	\describe{
	  \item{\code{Residence_Address_Longitude}}{Numeric vector, longitude coordinate of voter}
	  \item{\code{Residence_Address_Latitude}}{Numeric vector, latitude coordinate of voter}
	}
}
%\value{
%A dataset of 5,000 rows and two columns.
%}
\references{
King County, WA
}
\author{
Loren Collingwood <loren.collingwood@ucr.edu>
}
\examples{
data(king_dbox)
str(king_geo)
}