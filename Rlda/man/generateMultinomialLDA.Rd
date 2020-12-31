\name{generateMultinomialLDA}
\alias{generateMultinomialLDA}
\alias{generateMultinomialLDA.rlda}
\title{Simulates a Multinomial LDA.}
\description{
Simulates a Multinomial LDA.
}
\usage{
  \method{generateMultinomialLDA}{rlda}(seed0, community, variables,
  observations, totalElements, beta, gamma, ...)
}
\arguments{
  \item{seed0}{Initial seed to simulate a Multinomial LDA.}
  \item{community}{Total number of latent clusters. Must be greater than 2.}
  \item{variables}{Total number of variables. Must be greater than the number of communities.}
  \item{observations}{Total number of observations. Must be greater than 1.}
  \item{totalElements}{Total expected number of elements for each observation. Must be greater than 1.}
  \item{beta}{Vector of positive hyperparameters with dimension equal variables.}
  \item{gamma}{Scalar hyperparameters that must be positive.}
  \item{...}{ other arguments may be useful.}
}

\details{
Generates a list with the simulated Theta and Phi matrix of parameters, Z latent matrix of communities and and Data matrix for the Multinomial LDA.
}

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

\keyword{generateMultinomialLDA}
\keyword{generateBernoulliLDA}
\keyword{generateBinomialLDA}
\keyword{Rlda}

\seealso{\code{\link{generateBernoulliLDA}}, \code{\link{generateBinomialLDA}}}

\examples{
\dontrun{
#Generate fake data
res<- generateMultinomialLDA.rlda(seed0=5587, community=4,
         variables=100, observations=1000, totalElements=40,
         beta=rep(1,100), gamma=0.01)
#Show results
res
}
}
