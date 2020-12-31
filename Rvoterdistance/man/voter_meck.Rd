\name{voter_meck}
\alias{voter_meck}
\title{Dataset of registered voters, Mecklenburg County}
\description{Dataset of random registered voter locations in Mecklenburg County, North Carolina, as of 2016 general election.}
\usage{
data(meck_ev)
}
\format{A dataset with 4,552 rows and three columns:
	\describe{
		\item{\code{county}}{Character vector}
		\item{\code{long}}{Numeric vector, longitude coordinate}
		\item{\code{lat}}{Numeric vector, latitude coordinate}
	}
}
%\value{
%A dataset with 43 rows and five columns
%}
\references{
Mecklenburg County, North Carolina
}
\author{
Loren Collingwood <loren.collingwood@ucr.edu>, Hannah Walker <hlw25<at>georgetown.edu
}
\examples{
data(meck_ev) # Read in the stored RData file
str(voter_meck) # This is the actual dataset
}