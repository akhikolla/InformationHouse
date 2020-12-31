\name{generate.ESF}
\alias{generate.ESF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate community data under the standard neutral model of biodiversity, using the urn scheme as described in Etienne 2005
}
\description{
 This function generates community data under the standard neutral model of biodiversity, using the urn scheme as described in Etienne 2005
}
\usage{
generate.ESF(theta, I, J)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{Fundamental biodiversity number theta
}
  \item{I}{Fundamental dispersal number I
}
  \item{J}{total number of individuals in the local community
}
}
\value{Vector containing the unlabeled species abundances in the local community
}
\references{Etienne, R.S. (2005). A new sampling formula for neutral biodiversity. Ecology Letters, 8(3), 253-260.
}
\author{
	Thijs Janzen & Bart Haegeman
}
\examples{
	generate.ESF(theta=42,I=10,J=2000);
}
