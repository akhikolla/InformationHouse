\name{logLik}
\alias{logLik}
\alias{logLik.rlda}
\title{Provide the log-likelihood for the rlda object.}
\description{
Takes a \code{rlda} object produced and provides the log-likelihood.
}
\usage{
  \method{logLik}{rlda}(object, ...)
}
\arguments{
  \item{object}{a \code{rlda} object}
  \item{...}{ other arguments may be useful.}
}

\details{
Get the log-likelihood of the model.
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

\keyword{summary}
\keyword{LDA}
\keyword{Fuzzy}
\keyword{Cluster}

\seealso{\code{\link{rlda.binomial}}, \code{\link{rlda.bernoulli}},\code{\link{rlda.multinomial}} }

\examples{
\dontrun{
		library(Rlda)
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
		#log-likelihood
		logLik(res)
  }
}

