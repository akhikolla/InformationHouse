\name{summary}
\alias{summary}
\alias{summary.rlda}
\title{Provide summary information.}
\description{
Takes a \code{rlda} object produced by \code{rlda.binomial}, \code{rlda.bernoulli} or \code{rlda.multinomial} and produces a summary from it.
}
\usage{
  \method{summary}{rlda}(object, burnin=0.1, quantile = 0.95, silent=FALSE, ...)
}
\arguments{
  \item{object}{a \code{rlda} object as produced by \code{rlda.binomial}, \code{rlda.bernoulli} or \code{rlda.multinomial}.}
  \item{burnin}{a percentual of burn-in observations must be a number between 0 and 1. The default value is \code{burnin=0.1}}
 \item{quantile}{a quantile for the credible interval. The default value is \code{quantile=0.95}}
  \item{silent}{a logical value specifying if the number of clusters, variables and gibbs interations must be printed. The default value is \code{silent=FALSE}}
  \item{...}{ other arguments may be useful.}
}

\details{
Summarizes the Gibbs Samping results and arguments.
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
#Summary method
Theta<-summary(res)$Theta
}
}
