\name{EstimateIsing}
\alias{EstimateIsing}
\alias{EstimateIsingPL}
\alias{EstimateIsingUni}
\alias{EstimateIsingBi}
\alias{EstimateIsingLL}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
non-regularized estimation methods for the Ising Model
}
\description{
This function can be used for several non-regularized estimation methods of the Ising Model. See details.
}
\usage{
EstimateIsing(data, responses, beta = 1, method = c("pl", "uni",
                 "bi", "ll"), adj = matrix(1, ncol(data), ncol(data)),
                 ...)
EstimateIsingUni(data, responses, beta = 1, adj = matrix(1, ncol(data),
                 ncol(data)), ...)
EstimateIsingBi(data, responses, beta = 1, ...)
EstimateIsingPL(data, responses, beta = 1, ...)
EstimateIsingLL(data, responses, beta = 1, adj = matrix(1, ncol(data),
                 ncol(data)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
 Data frame with binary responses to estimate the Ising model over
}
  \item{responses}{
 Vector of length two indicating the response coding (usually \code{c(0L, 1L)} pr \code{c(-1L, 1L)})
}
  \item{beta}{
 Inverse temperature parameter
}
  \item{method}{
  The method to be used. \code{pl} uses pseudolikelihood estimation, \code{uni} sequential univariate regressions, \code{bi} bivariate regressions and \code{ll} estimates the Ising model as a loglinear model.
}
\item{adj}{
Adjacency matrix of the Ising model.
}
  \item{\dots}{
 Arguments sent to estimator functions
}
}
\details{
The following algorithms can be used (see Epskamp, Maris, Waldorp, Borsboom; in press).
\describe{
\item{\code{pl}}{Estimates the Ising model by maximizing the pseudolikelihood (Besag, 1975).}
\item{\code{uni}}{Estimates the Ising model by computing univariate logistic regressions of each node on all other nodes. This leads to a single estimate for each threshold and two estimates for each network parameter. The two estimates are averaged to produce the final network. Uses \code{\link{glm}}.}
\item{\code{bi}}{Estimates the Ising model using multinomial logistic regression of each pair of nodes on all other nodes. This leads to a single estimate of each network parameter and $p$ estimates of each threshold parameter. Uses \code{\link[nnet]{multinom}}.}
\item{\code{ll}}{Estimates the Ising model by phrasing it as a loglinear model with at most pairwise interactions. Uses \code{\link{loglin}}.}
}
}
\value{
A list containing the estimation results:
\item{graph}{The estimated network}
\item{thresholds}{The estimated thresholds}
\item{results}{The results object used in the analysis}
%% ...
}
\references{
Epskamp, S., Maris, G., Waldorp, L. J., and Borsboom, D. (in press). Network Psychometrics. To appear in: Irwing, P., Hughes, D., and Booth, T. (Eds.), Handbook of Psychometrics. New York: Wiley.

Besag, J. (1975), Statistical analysis of non-lattice data. The statistician, 24, 179-195.
}
\author{
Sacha Epskamp (mail@sachaepskamp.com)
}
\examples{
# Input:
N <- 5 # Number of nodes
nSample <- 500 # Number of samples

# Ising parameters:
Graph <- matrix(sample(0:1,N^2,TRUE,prob = c(0.7, 0.3)),N,N) * rnorm(N^2)
Graph <- Graph + t(Graph)
diag(Graph) <- 0
Thresholds <- rep(0,N)
Beta <- 1

# Response options (0,1 or -1,1):
Resp <- c(0L,1L)
Data <- IsingSampler(nSample,Graph, Thresholds)

# Pseudolikelihood:
resPL <- EstimateIsing(Data, method = "pl")
cor(Graph[upper.tri(Graph)], resPL$graph[upper.tri(resPL$graph)])

# Univariate logistic regressions:
resUni <- EstimateIsing(Data, method = "uni")
cor(Graph[upper.tri(Graph)], resUni$graph[upper.tri(resUni$graph)])

# bivariate logistic regressions:
resBi <- EstimateIsing(Data, method = "bi")
cor(Graph[upper.tri(Graph)], resBi$graph[upper.tri(resBi$graph)])

# Loglinear model:
resLL <- EstimateIsing(Data, method = "ll")
cor(Graph[upper.tri(Graph)], resLL$graph[upper.tri(resLL$graph)])
}
