\name{rlda.binomial}
\alias{rlda.binomial}
\title{LDA with binomial entry and Stick-Breaking prior.}
\description{
  This method implements the Latent Dirichlet Allocation with
  Stick-Breaking prior for binomial data.
  \code{rlda.binomial} works with frequency data.frame and also a
  population data.frame.
}
\usage{
  rlda.binomial(data, pop, n_community, alpha0, alpha1, gamma,
  n_gibbs, ll_prior = TRUE, display_progress = TRUE)
}
\arguments{
  \item{data}{A abundance data.frame where each row is a sampling unit
    (i.e. Plots, Locations, Time, etc.) and each column is a categorical
    type of element (i.e. Species, Firms, Issues, etc.).}
  \item{pop}{A population data.frame where each row is a sampling unit
  (i.e. Plots, Locations, Time, etc.) and each column is a categorical
  type of element (i.e. Species, Firms, Issues, etc.). The elements inside
  this data.frame must all be greater than the elements inside the \code{data}
  data.frame.}
  \item{n_community}{Total number of communities to return. It must be less than
  the total number of columns inside the \code{data} and \code{pop} data.frame.}
  \item{alpha0}{Hyperparameter associated with the Beta prior Beta(alpha0, alpha1).}
  \item{alpha1}{Hyperparameter associated with the Beta prior Beta(alpha0, alpha1).}
  \item{gamma}{Hyperparameter associated with the Stick-Breaking prior.}
  \item{n_gibbs}{Total number of Gibbs Samples.}
  \item{ll_prior}{boolean scalar, \code{TRUE} if the log-likelihood
  must be computed using also the priors or \code{FALSE} otherwise.}
  \item{display_progress}{boolean scalar, \code{TRUE} if the
  Progress Bar must be showed and \code{FALSE} otherwise.}
      }
    \value{
      A R List with three elements:
  \item{Theta}{The individual probability for each observation
  (ex: location) belong in each cluster (ex: community). It is a matrix
  with dimension equal \code{n_gibbs} by \code{nrow(data) * n_community}}
  \item{Phi}{The individual probability for each variable
  (ex: Specie) belong in each cluster (ex: community). It is a matrix
  with dimension equal \code{n_gibbs} by \code{ncol(data) * n_community}}
  \item{LogLikelihood}{The vector of Log-Likelihoods compute for
  each Gibbs Sample.}
}
\details{
  \code{rlda.binomial} uses a modified Latent Dirichlet Allocation method
  to construct Mixed-Membership Clusters using Bayesian Inference.
  The \code{data} must be a non-empty data.frame with the frequencies for each variable
  (column) in each observation (row). The \code{pop} must be a non-empty data.frame with
  the frequencies for each variable (column) in each observation (row) greater than the
  entries inside \code{data} data.frame.
}
\note{The \code{Theta} and \code{Phi} matrix can be obtained for the i-th gibbs
sampling using \code{matrix(Theta[i,], nrow = nrow(data), ncol = n_community)} and
\code{matrix(Phi[i,], nrow = n_community, ncol = ncol(data))}, respectively.}
\references{
  \itemize{
    \item
     Blei, David M., Andrew Y. Ng, and Michael I. Jordan.\cr
     "Latent dirichlet allocation." \emph{Journal of machine Learning research}
      3.Jan (2003): 993-1022.
      \url{http://www.jmlr.org/papers/volume3/blei03a/blei03a.pdf}
   \item
      Valle, Denis, et al.\cr
      "Decomposing biodiversity data using the Latent Dirichlet
      Allocation model, a probabilistic multivariate statistical
      method." \emph{Ecology letters} 17.12 (2014): 1591-1601.
  }
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
      \url{http://daijiang.name}
  }
}
\seealso{\code{\link{rlda.multinomial}}, \code{\link{rlda.bernoulli}} }

\keyword{LDA}
\keyword{Mixed-Membership}
\keyword{Cluster}
\keyword{Binomial}

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
	}
}

