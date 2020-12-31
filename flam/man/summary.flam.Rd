\name{summary.flam}
\alias{summary.flam}
\alias{plot.flamSparsity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes a Call to \code{flam}
}
\description{
This function summarizes a call to \code{\link{flam}}, as well as the sparsity pattern of the resulting feature estimates for a single or all fits.
}
\usage{
\method{summary}{flam}(object, index = NULL, \dots)
\method{plot}{flamSparsity}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "flam".
}
  \item{index}{
the index for the fit of interest to be summarized. Note that \code{index} of i corresponds to the model with tuning parameters \code{object$all.alpha[i]} and \code{object$all.lambda[i]}. If \code{index} is not specified, information summarizing all fits is given.
}
  \item{x}{
an object of class 'flamSparsity', which is silently returned by \code{summary.flam}. 
}
  \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}

\value{
	If \code{index} is not specified, \code{summary.flam} silently returns the sparsity matrix and tuning parameters in an object of class 'flamSparsity'. This is used when \code{plot(summary(object))} is called.
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
data <- sim.data(n = 50, scenario = 1, zerof = 10, noise = 1)
#fit model for a range of tuning parameters
flam.out <- flam(x = data$x, y = data$y, alpha.seq = c(0.8, 0.9, 1))

#summarize all of the fits (in particular, the percent sparsity achieved)
#note: percent sparsity is the percentage of features estimated to have 
#no relationship with outcome
summary(flam.out)
#the percent sparsity over the range of tuning parameters can also 
#be displayed in a plot
plot(summary(flam.out))

#we can also summarize the fit with a specific index, e.g. 25
#that is, lambda and alpha of
flam.out$all.lambda[25]; flam.out$all.alpha[25]
summary(flam.out, index = 25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
