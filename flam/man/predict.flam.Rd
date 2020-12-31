\name{predict.flam}
\alias{predict.flam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicts Observations for a New Covariate Matrix and Fit from \code{\link{flam}}
}
\description{
This function makes predictions from a specified covariate matrix for a fit of the class "flam" with user-specified tuning parameters.
}
\usage{
\method{predict}{flam}(object, new.x, lambda, alpha, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of the class "flam".
}
  \item{new.x}{
the covariate matrix for which to make predictions - the number of columns should match that of \code{object$x}.
}
  \item{lambda}{
the desired value for the tuning parameter lambda. This does not need to be a value in \code{object$all.lambda}.
}
  \item{alpha}{
the desired value for the tuning parameter alpha. This does not need to be a value in \code{object$all.alpha}.
}
  \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}
\details{
It is likely that \code{new.x[,i]} contains values not contained in \code{object$x[,i]}. Predictions for that particular case are taken to be a linear interpolation of the nearest neighboring values in \code{object$x[,i]}, i.e., the closest smaller value and the closest larger value.
}
\value{
A vector containing the fitted y values for \code{new.x}.
}
\references{
Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\author{
Ashley Petersen
}

\seealso{
\code{\link{flam}}
}
\examples{
#See ?'flam-package' for a full example of how to use this package

#generate data
set.seed(1)
data <- sim.data(n = 100, scenario = 1, zerof = 0, noise = 1)

#fit model for a range of tuning parameters
flam.out <- flam(x = data$x, y = data$y)

#we can make predictions for a covariate matrix with new observations
#choose desired alpha and lambda
alpha <- flam.out$all.alpha[15]; lambda <- flam.out$all.lambda[15]
#new.x with 20 observations and the same number of features as flam.out$x
new.data <- sim.data(n = 20, scenario = 1, zerof = 0, noise = 1)
new.x <- new.data$x
#make predictions
y.hat <- predict(flam.out, new.x = new.x, lambda = lambda, alpha = alpha)
#which can be compared to the true y
plot(new.data$y, y.hat, xlab="y", ylab=expression(hat(y)))
abline(0,1,lty=2)

#can also make predictions for any alpha and lambda:
predict(flam.out, new.x = new.x, lambda = 2, alpha = 0.9)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
