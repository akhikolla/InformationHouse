\name{plot}
\alias{plot.rlda}
\alias{plot,Rlda}
\title{plot method for rlda object}


\description{Plot a rlda object.
The \code{plot} function returns three plots based on Theta matrix, Phi matrix and log-likelihood. }

\usage{
  \S3method{plot}{rlda}(x, burnin=0.1, maxCluster=NA, ...)
}

\arguments{

  \item{x}{a \code{rlda} object created by \code{rlda.binomial}, \code{rlda.bernoulli} or \code{rlda.multinomial} function.}
  \item{burnin}{a percentual of burn-in observations must be a number between 0 and 1. The default value is \code{burnin=0.1}.}
  \item{maxCluster}{The maximum number of cluster to be shown. The default value is \code{maxCluster=NA} which represents all clusters must be shown.}
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

\keyword{plot}
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
#Plot the results
plot(res)
}
}
