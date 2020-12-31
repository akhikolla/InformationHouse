\name{summary.flamCV}
\alias{summary.flamCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes a Call to \code{flamCV}
}
\description{
This function summarizes a call to code{\link{flamCV}} and identifies the tuning parameter chosen by cross-validation.
}
\usage{
\method{summary}{flamCV}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "flamCV".
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
flamCV.out <- flamCV(x = data$x, y = data$y, n.fold = 2)

#we can summarize the cross-validation function call
summary(flamCV.out)
#lambda chosen by cross-validation is also available from
flamCV.out$lambda.cv
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
