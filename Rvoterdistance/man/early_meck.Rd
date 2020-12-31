
\name{early_meck}
\alias{early_meck}
\title{Dataset of early vote locations}
\description{Dataset of early vote locations in Mecklenburg County, North Carolina, as of 2016 general election.}
\usage{
data(meck_ev)
}
\format{A dataset with 21 rows and five columns:
	\describe{
		\item{\code{match_addr}}{Character vector}
		\item{\code{county}}{Character vector}
		\item{\code{office}}{Character vector}
		\item{\code{long}}{Numeric vector, longitude coordinate}
		\item{\code{lat}}{Numeric vector, latitude coordinate}
	}
}
\references{
Mecklenburg County, Loren Collingwood
}
\author{
Loren Collingwood <loren.collingwood@ucr.edu>, Hannah Walker <hlw25<at>georgetown.edu
}
\examples{
data(meck_ev)
str(early_meck)
}