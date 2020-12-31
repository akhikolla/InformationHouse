\name{rlda.multinomial}
\alias{rlda.multinomial}
\title{LDA with multinomial entry and Stick-Breaking prior.}
\description{
  This method implements the Latent Dirichlet Allocation with
  Stick-Breaking prior for multinomial data.
  \code{rlda.multinomial} works with frequency data.frame.
}
\usage{
  rlda.multinomial(data, n_community, beta, gamma,
  n_gibbs, ll_prior = TRUE, display_progress = TRUE)
}
\arguments{
  \item{data}{A abundance data.frame where each row is a sampling unit
    (i.e. Plots, Locations, Time, etc.) and each column is a categorical
    type of element (i.e. Species, Firms, Issues, etc.).}
  \item{n_community}{Total number of communities to return. It must be less than
  the total number of columns inside the \code{data} data.frame.}
  \item{beta}{Hyperparameter associated with the Dirichlet \code{Phi}
  matrix.}
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
  \code{rlda.multinomial} uses a modified Latent Dirichlet Allocation method
  to construct Mixed-Membership Clusters using Bayesian Inference.
  The \code{data} must be a non-empty data.frame with the frequencies for each variable
  (column) in each observation (row).
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
  }
}
\seealso{\code{\link{rlda.binomial}}, \code{\link{rlda.bernoulli}} }

\keyword{LDA}
\keyword{Mixed-Membership}
\keyword{Cluster}
\keyword{Multinomial}

\examples{
	\dontrun{
		# Invoke the library
		library(Rlda)
		# Read the Complaints data
		data(complaints)

		# Create the abundance matrix
		library(reshape2)
		mat1 <- dcast(complaints[, c("Company","Issue")],
					Company ~ Issue, fun.aggregate = length,
					value.var = "Issue")
		# Create the rowname
		rownames(mat1) <- mat1[, 1]
		# Remove the ID variable
		mat1 <- mat1[, -1]

		# Set seed
		set.seed(9292)
		# Hyperparameters for each prior distribution
		beta <- rep(1,ncol(mat1))
		gamma <- 0.01

		#Execute the LDA for the Multinomial entry
		res <- rlda.multinomial(data = mat1, n_community = 30,
		beta = beta, gamma = gamma, n_gibbs = 1000,
		ll_prior = TRUE, display_progress = TRUE)
	}
}

