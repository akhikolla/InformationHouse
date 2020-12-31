\name{expected.SAD}
\alias{expected.SAD}
\title{
  Calculate the expected species abundance distribution of the standard neutral model, given theta, m and J
}
\description{
  This function calculates the expected species abundance distribution of the standard neutral model given theta, m and J, sensu equation 6 from Etienne and Alonso (2005). 
}
\usage{
expected.SAD(theta, m, J)
}

\arguments{
  \item{theta}{Fundamental biodiversity number theta

}
  \item{m}{migration parameter

}
  \item{J}{Total number of individuals in the local community

}
}

\value{A vector containing the abundances binned into log2 bins (sensu Preston).
}
\references{
	Etienne, R.S., & Alonso, D. (2005). A dispersal-limited sampling theory for species and alleles. Ecology Letters, 8(100), 1147-1156.
}
\author{
	Thijs Janzen & Bart Haegeman
}

\examples{
SAD <- expected.SAD(theta=42, m = 0.1, J = 200);
barplot(SAD,names.arg=0:(length(SAD)-1),
	xlab="Number of individuals (log2)",
	ylab="Number of Species" )
}