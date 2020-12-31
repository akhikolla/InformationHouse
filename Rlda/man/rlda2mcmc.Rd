\name{rlda2mcmc}
\alias{rlda2mcmc}
\alias{rlda2mcmc.rlda}
\title{Conversion between rlda to mcmc object from coda package}
\description{Convert rlda object to mcmc object from coda package
The \code{rlda2mcmc} function create a list of two elements incluing the Theta matrix and Phi matrix in a mcmc class.
}
\usage{
  \method{rlda2mcmc}{rlda}(object, ...)
}
\arguments{
  \item{object}{a \code{rlda} object created by \code{rlda.binomial}, \code{rlda.bernoulli} or \code{rlda.multinomial} function.}
  \item{...}{ other arguments may be useful.}
}

\seealso{\code{\link{rlda.binomial}}, \code{\link{rlda.bernoulli}},\code{\link{rlda.multinomial}} }

\author{
  \itemize{
    \item
      Pedro Albuquerque.\cr
      \email{pedroa@unb.br}\cr
      \url{http://pedrounb.blogspot.com/}
   \item
      Denis Valle.\cr
      \email{drvalle@ufl.edu}\cr
      \url{http://denisvalle.weebly.com/}
   \item
      Daijiang Li.\cr
      \email{daijianglee@gmail.com}\cr
  }
}

\keyword{rlda2mcmc}
\keyword{LDA}
\keyword{Fuzzy}
\keyword{Cluster}

\examples{
\dontrun{
#Load data
data(presence)
#Set seed
set.seed(9842)
#Hyperparameters for each prior distribution
gamma <-0.01
alpha0<-0.01
alpha1<-0.01
#Execute the LDA for the Binomial entry
res<-rlda.bernoulli(presence, 10, alpha0, alpha1, gamma,
                    5000, TRUE, FALSE)
#rlda2mcmc the results
listParams <- rlda2mcmc(res)
}
}
