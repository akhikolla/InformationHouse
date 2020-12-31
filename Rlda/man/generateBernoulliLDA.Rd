\name{generateBernoulliLDA}
\alias{generateBernoulliLDA}
\alias{generateBernoulliLDA.rlda}
\title{Simulates a Bernoulli LDA.}
\description{
Simulates a Bernoulli LDA.
}
\usage{
  \method{generateBernoulliLDA}{rlda}(seed0, community, variables,
  observations, alpha0, alpha1, gamma, ...)
}
\arguments{
  \item{seed0}{Initial seed to simulate a Bernoulli LDA.}
  \item{community}{Total number of latent clusters. Must be greater than 2.}
  \item{variables}{Total number of variables. Must be greater than the number of communities.}
  \item{observations}{Total number of observations. Must be greater than 1.}
  \item{alpha0}{Scalar hyperparameters that must be positive.}
  \item{alpha1}{Scalar hyperparameters that must be positive.}
  \item{gamma}{Scalar hyperparameters that must be positive.}
  \item{...}{ other arguments may be useful.}
}
\details{
Generates a list with the simulated Theta and Phi matrix of parameters,
Z latent matrix of communities and and Data matrix for the Bernoulli LDA.
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

\seealso{\code{\link{generateMultinomialLDA}}, \code{\link{generateBinomialLDA}}}

\examples{
\dontrun{
#Generate fake data
res<- generateBernoulliLDA.rlda(seed0=9292, community=3,
                                variables=100, observations=1000,
                                alpha0=0.01, alpha1=0.01, gamma=0.01)
#Show results
res
}
}
