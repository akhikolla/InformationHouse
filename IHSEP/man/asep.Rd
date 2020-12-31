\name{asep}
\alias{asep}
\docType{data}
\title{
An IHSEP data set
}
\description{
A simulated data set from the inhomegeneous self-exciting process model
with baseline intensity, or immigration rate,
\eqn{\nu(t)=200(2+\cos(2\pi t))} and excitation function, or fertility
\eqn{g(t)=8\exp(-16 t)}. 
}
\usage{data(asep)}
\format{
  The format is a list of the arrival/birth times of individuals of all
  generations, in the order of generation 0 individuals, generation 1
  individuals, etc. 
}
\details{
  Times of arrivals/births of the same generation is listed
  together in ascending order. Number of generations is given by the
  length of the object. 
}
\source{
  Simulated by a call to the function \code{simHawkes1}.
}
\examples{
data(asep)
## number of generations
length(asep)
## jump times of the observed self-exciting process
sort(unlist(asep))
}
\keyword{datasets}
