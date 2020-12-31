\name{flam}
\alias{flam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the Fused Lasso Additive Model for a Sequence of Tuning Parameters
}
\description{
Fit an additive model where each component is estimated to piecewise constant with a small number of adaptively-chosen knots. The model is fit for a sequence of tuning parameters. In particular, this function implements the "fused lasso additive model", as proposed in Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\usage{
flam(x, y, lambda.min.ratio = 0.01, n.lambda = 50, lambda.seq = NULL,
alpha.seq = 1, family = "gaussian", method = "BCD", tolerance = 10e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
n x p covariate matrix. May have p > n.
}
  \item{y}{
n-vector containing the outcomes for the n observations in \code{x}.
}
  \item{lambda.min.ratio}{
smallest value for \code{lambda.seq}, as a fraction of the maximum lambda value, which is the data-derived smallest value for which all estimated functions are zero. The default is 0.01.
}
  \item{n.lambda}{
the number of lambda values to consider - the default is 50.
}
  \item{lambda.seq}{
a user-supplied sequence of positive lambda values to consider. The typical usage is to calculate \code{lambda.seq} using \code{lambda.min.ratio} and \code{n.lambda}, but providing \code{lambda.seq} overrides this. If provided, \code{lambda.seq} should be a decreasing sequence of values, since \code{flam} relies on warm starts for speed. Thus fitting the model for a whole sequence of lambda values is often faster than fitting for a single lambda value. Note that the model is fit for all combinations of \code{alpha.seq} and \code{lambda.seq}, so all values of \code{lambda.seq} provided should be unique.
}
  \item{alpha.seq}{
the value(s) of alpha to consider - default is 1. Values must be in [0,1] with values near 0 prioritizing sparsity of functions and values near 1 prioritizing limiting the number of knots. Empirical evidence suggests using alpha of 1 when p < n and alpha of 0.75 when p > n. Note that the model is fit for all combinations of \code{alpha.seq} and \code{lambda.seq}, so all values of \code{alpha.seq} provided should be unique.
}
  \item{family}{
specifies the loss function to use. Currently supports squared error loss (default; \code{family="gaussian"}) and logistic loss (\code{family="binomial"}).
}
  \item{method}{
specifies the optimization algorithm to use. Options are block-coordinate descent (default; \code{method="BCD"}), generalized gradient descent (\code{method="GGD"}), or generalized gradient descent with backtracking (\code{method="GGD.backtrack"}). This argument is ignored if \code{family="binomial"}.
}
  \item{tolerance}{
specifies the convergence criterion for the objective (default is 10e-6).
}
}

\value{
An object with S3 class "flam".
\item{all.alpha}{vector of alpha values considered. This will be m times longer than the user-specified \code{alpha.seq} where m is the length of the user-specified \code{lambda.seq}.}
\item{all.lambda}{vector of lambda values considered. This will be q times longer than the user-specified \code{lambda.seq} where q is the length of the user-specified \code{alpha.seq}.}
\item{theta.hat.list}{list of estimated theta matrices of dimension n x p. Note that the predicted values \code{y.hat.mat[i,] = g(}\code{beta0.hat.vec[i] + }\code{rowSums(theta.hat.list[[i]]))} where \code{g} is the link function (identity if \code{family="gaussian"} and expit if \code{family="binomial"}).}
\item{f.hat.list}{list of estimated function matrices of dimension n x p. Note that \code{f.hat.list[[i]]} is \code{theta.hat.list[[i]]} with the elements of each column ordered in terms of increasing \code{x[,i]}.}
\item{beta0.hat.vec}{vector of estimated intercepts with \code{beta0.hat.vec[i]} being the intercept for the model with tuning parameters \code{all.alpha[i]} and \code{all.lambda[i]}.}
\item{y.hat.mat}{matrix with \code{y.hat.mat[i,]} containing fitted y values for the tuning parameters \code{all.alpha[i]} and \code{all.lambda[i]}.}
\item{non.sparse.list}{list with \code{non.sparse.list[[i]]} containing the indices for the predictors with non-sparse fits for the tuning parameters \code{all.alpha[i]} and \code{all.lambda[i]}.}
\item{num.non.sparse}{vector with \code{num.non.sparse[i]} indicating the number of non-sparse predictor fits for the tuning parameters \code{all.alpha[i]} and \code{all.lambda[i]}.}
\item{y}{as specified by user.}
\item{x}{as specified by user.}
\item{family}{as specified by user (or default).}
\item{method}{as specified by user (or default).}
\item{tolerance}{as specified by user (or default).}
\item{call}{the matched call.}
}
\references{
Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\author{
Ashley Petersen
}

\seealso{
\code{\link{predict.flam}}, \code{\link{plot.flam}}, \code{\link{summary.flam}}
}
\examples{
#See ?'flam-package' for a full example of how to use this package

#generate data
set.seed(1)
data <- sim.data(n = 50, scenario = 1, zerof = 10, noise = 1)

#fit model for a range of lambda chosen by default and alpha's of 0.75 and 1
flam.out <- flam(x = data$x, y = data$y, alpha.seq = c(0.75, 1))
#or specify desired lambda sequence (often equally spaced on log scale)
#should be a decreasing sequence of several values for computational speed
user.lambda.seq <- exp(seq(log(50), log(1), len=40))
flam.out2 <- flam(x = data$x, y = data$y, lambda.seq = user.lambda.seq)

\dontrun{
#alternatively, generate data for logistic FLAM model
data2 <- sim.data(n = 50, scenario = 1, zerof = 10, family = "binomial")
#fit the FLAM model using logistic loss
flam.logistic.out <- flam(x = data2$x, y = data2$y, family = "binomial")
}

#'flam' returns an object of the class 'flam'
#see ?'flam-package' for an example using S3 methods for 'flam' objects

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
