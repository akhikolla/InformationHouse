
\name{dbox}
\alias{dbox}
\title{Dataset of drop box locations}
\description{Dataset of drop box locations in King County, Washington, as of 2016 general election.}
\usage{
data(king_dbox)
}
\format{A dataset with 43 rows and five columns:
	\describe{
		\item{\code{location_name}}{Character vector}
		\item{\code{address_city}}{Character vector}
		\item{\code{state}}{Character vector}
		\item{\code{lat}}{Numeric vector, latitude coordinate}
		\item{\code{long}}{Numeric vector, longitude coordinate}
	}
}
%\value{
%A dataset with 43 rows and five columns
%}
\references{
King County, Washington
}
\author{
Loren Collingwood <loren.collingwood@ucr.edu>
}
\examples{
data(king_dbox)
str(dbox)
}