\name{plot.flamCV}
\alias{plot.flamCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots Cross-Validation Curve for Object of Class "flamCV"
}
\description{
This function plots the cross-validation curve for a series of models fit using \code{\link{flamCV}}. The cross-validation error with +/-1 standard error is plotted for each value of lambda considered in the call to \code{flamCV} with a dotted vertical line indicating the chosen lambda.
}
\usage{
\method{plot}{flamCV}(x, showSE = T, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class "flamCV".
}
  \item{showSE}{
a logical (\code{TRUE} or \code{FALSE}) for whether the standard errors of the curve should be plotted.
}
  \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}
\references{
Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\author{
Ashley Petersen
}

\seealso{
\code{\link{flamCV}}
}
\examples{
#See ?'flam-package' for a full example of how to use this package

#generate data
set.seed(1)
data <- sim.data(n = 50, scenario = 1, zerof = 0, noise = 1)

#fit model and select tuning parameters using 2-fold cross-validation
#note: use larger 'n.fold' (e.g., 10) in practice
flamCV.out <- flamCV(x = data$x, y = data$y, within1SE = TRUE, n.fold = 2)

#lambdas chosen is
flamCV.out$lambda.cv

#we can now plot the cross-validation error curve with standard errors
#vertical dotted line at lambda chosen by cross-validation
plot(flamCV.out)
#or without standard errors
plot(flamCV.out, showSE = FALSE)

\dontrun{
#can choose lambda to be value with minimum CV error
#instead of lambda with CV error within 1 standard error of the minimum
flamCV.out2 <- flamCV(x = data$x, y = data$y, within1SE = FALSE, n.fold = 2)

#contrast to chosen lambda for minimum cross-validation error
#it's a less-regularized model (i.e., lambda is smaller)
plot(flamCV.out2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
