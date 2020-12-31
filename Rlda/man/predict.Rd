\name{predict}
\alias{predict}
\alias{predict.rlda}
\title{Provide predictions to the Binomial entry.}
\description{
Takes a \code{rlda} object produced by \code{rlda.binomial} and produces a prediction from it.
}
\usage{
  \method{predict}{rlda}(object, data, nclus=NA, burnin=0.1, places.round=0, ...)
}
\arguments{
  \item{object}{a \code{rlda} object as produced by \code{rlda.binomial}}
  \item{data}{Dataset used to make the predictions. Must have the same number of columns as the dataset used in the \code{rlda.binomial}.}
  \item{nclus}{Number of clusters to be used in the prediction. The default value is \code{nclus=NA}}
  \item{burnin}{a percentual of burn-in observations must be a number between 0 and 1. The default value is \code{burnin=0.1}}
  \item{places.round}{Number decimal places tob rounded. The default value is \code{places.round=0}}
  \item{...}{ other arguments may be useful.}
}

\details{
Predicts the Gibbs Samping results and arguments.
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
		# Read the SP500 data
		data(sp500)
		# Create size
		spSize <- as.data.frame(matrix(100,
							  ncol = ncol(sp500),
							  nrow = nrow(sp500)))
		# Set seed
		set.seed(5874)
		# Hyperparameters for each prior distribution
		gamma  <- 0.01
		alpha0 <- 0.01
		alpha1 <- 0.01
		# Execute the LDA for the Binomial entry
		res <- rlda.binomial(data = sp500, pop = spSize, n_community = 10,
		alpha0 = alpha0, alpha1 = alpha1, gamma = gamma,
		n_gibbs = 500, ll_prior = TRUE, display_progress = TRUE)
		#Predict
		pred<- predict(res, sp500, nclus=3)
  }
}
