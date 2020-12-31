\name{flamCV}
\alias{flamCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the Fused Lasso Additive Model and Do Tuning Parameter Selection using K-Fold Cross-Validation
}
\description{
Fit an additive model where each component is estimated to piecewise constant with a small number of adaptively-chosen knots. Tuning parameter selection is done using K-fold cross-validation. In particular, this function implements the "fused lasso additive model", as proposed in Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\usage{
flamCV(x, y, lambda.min.ratio = 0.01, n.lambda = 50, lambda.seq = NULL,
alpha = 1, family = "gaussian", method = "BCD", fold = NULL,
n.fold = NULL, seed = NULL, within1SE = T, tolerance = 10e-6)
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
a user-supplied sequence of positive lambda values to consider. The typical usage is to calculate \code{lambda.seq} using \code{lambda.min.ratio} and \code{n.lambda}, but providing \code{lambda.seq} overrides this. If provided, \code{lambda.seq} should be a decreasing sequence of values, since \code{flamCV} relies on warm starts for speed. Thus fitting the model for a whole sequence of lambda values is often faster than fitting for a single lambda value.
}
  \item{alpha}{
the value of the tuning parameter alpha to consider - default is 1. Value must be in [0,1] with values near 0 prioritizing sparsity of functions and values near 1 prioritizing limiting the number of knots. Empirical evidence suggests using alpha of 1 when p < n and alpha of 0.75 when p > n.
}
  \item{family}{
specifies the loss function to use. Currently supports squared error loss (default; \code{family="gaussian"}) and logistic loss (\code{family="binomial"}).
}
  \item{method}{
specifies the optimization algorithm to use. Options are block-coordinate descent (default; \code{method="BCD"}), generalized gradient descent (\code{method="GGD"}), or generalized gradient descent with backtracking (\code{method="GGD.backtrack"}).  This argument is ignored if \code{family="binomial"}.
}
  \item{fold}{
user-supplied fold numbers for cross-validation. If supplied, \code{fold} should be an n-vector with entries in 1,...,K when doing K-fold cross-validation. The default is to choose \code{fold} using \code{n.fold}.
}
  \item{n.fold}{
the number of folds, K, to use for the K-fold cross-validation selection of tuning parameters. The default is 10 - specification of \code{fold} overrides use of \code{n.fold}.
}
  \item{seed}{
an optional number used with \code{set.seed()} at the beginning of the function. This is only relevant if \code{fold} is not specified by the user.
}
  \item{within1SE}{
logical (\code{TRUE} or \code{FALSE}) for how cross-validated tuning parameters should be chosen. If \code{within1SE=TRUE}, lambda is chosen to be the value corresponding to the most sparse model with cross-validation error within one standard error of the minimum cross-validation error. If \code{within1SE=FALSE}, lambda is chosen to be the value corresponding to the minimum cross-validation error.
}
  \item{tolerance}{
specifies the convergence criterion for the objective (default is 10e-6).
}
}
\details{
Note that \code{flamCV} does not cross-validate over \code{alpha} - just a single value should be provided. However, if the user would like to cross-validate over \code{alpha}, then \code{flamCV} should be called multiple times for different values of \code{alpha} and the same \code{seed}. This ensures that the cross-validation folds (\code{fold}) remain the same for the different values of \code{alpha}. See the example below for details.
}
\value{
An object with S3 class "flamCV".
\item{mean.cv.error}{m-vector containing cross-validation error where m is the length of \code{lambda.seq}. Note that \code{mean.cv.error[i]} contains the cross-validation error for tuning parameters \code{alpha} and \code{flam.out$all.lambda[i]}.}
\item{se.cv.error}{m-vector containing cross-validation standard error where m is the length of \code{lambda.seq}. Note that \code{se.cv.error[i]} contains the standard error of the cross-validation error for tuning parameters \code{alpha} and \code{flam.out$all.lambda[i]}.}
\item{lambda.cv}{optimal lambda value chosen by cross-validation.}
\item{alpha}{as specified by user (or default).}
\item{index.cv}{index of the model corresponding to 'lambda.cv'.}
\item{flam.out}{object of class 'flam' returned by \code{flam}.}
\item{fold}{as specified by user (or default).}
\item{n.folds}{as specified by user (or default).}
\item{within1SE}{as specified by user (or default).}
\item{tolerance}{as specified by user (or default).}
\item{call}{matched call.}
}
\references{
Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\author{
Ashley Petersen
}

\seealso{
\code{\link{flam}}, \code{\link{plot.flamCV}}, \code{\link{summary.flamCV}}
}
\examples{
#See ?'flam-package' for a full example of how to use this package

#generate data
set.seed(1)
data <- sim.data(n = 50, scenario = 1, zerof = 10, noise = 1)

#fit model for a range of lambda chosen by default
#pick lambda using 2-fold cross-validation
#note: use larger 'n.fold' (e.g., 10) in practice
flamCV.out <- flamCV(x = data$x, y = data$y, alpha = 0.75, n.fold = 2)

\dontrun{
#note that cross-validation is only done to choose lambda for specified alpha
#to cross-validate over alpha also, call 'flamCV' for several alpha and set seed
#note: use larger 'n.fold' (e.g., 10) in practice
flamCV.out1 <- flamCV(x = data$x, y = data$y, alpha = 0.65, seed = 100, 
	within1SE = FALSE, n.fold = 2)
flamCV.out2 <- flamCV(x = data$x, y = data$y, alpha = 0.75, seed = 100, 
	within1SE = FALSE, n.fold = 2)
flamCV.out3 <- flamCV(x = data$x, y = data$y, alpha = 0.85, seed = 100, 
	within1SE = FALSE, n.fold = 2)
#this ensures that the folds used are the same
flamCV.out1$fold; flamCV.out2$fold; flamCV.out3$fold
#compare the CV error for the optimum lambda of each alpha to choose alpha
CVerrors <- c(flamCV.out1$mean.cv.error[flamCV.out1$index.cv], 
	flamCV.out2$mean.cv.error[flamCV.out2$index.cv], 
	flamCV.out3$mean.cv.error[flamCV.out3$index.cv])
best.alpha <- c(flamCV.out1$alpha, flamCV.out2$alpha, 
	flamCV.out3$alpha)[which(CVerrors==min(CVerrors))]

#also can generate data for logistic FLAM model
data2 <- sim.data(n = 50, scenario = 1, zerof = 10, family = "binomial")
#fit the FLAM model with cross-validation using logistic loss
#note: use larger 'n.fold' (e.g., 10) in practice
flamCV.logistic.out <- flamCV(x = data2$x, y = data2$y, family = "binomial",
	n.fold = 2)
}

#'flamCV' returns an object of the class 'flamCV' that includes an object
#of class 'flam' (flam.out); see ?'flam-package' for an example using S3
#methods for the classes of 'flam' and 'flamCV'
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
