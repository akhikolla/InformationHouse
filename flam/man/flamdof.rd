\name{flamDOF}
\alias{flamDOF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Degrees of Freedom for Fused Lasso Additive Model
}
\description{
This function calculates the degrees of freedom for a fused lasso additive model fit using \code{\link{flam}}.
}
\usage{
flamDOF(object, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of the class "flam".
}
  \item{index}{
the index for the model of interest. Note that \code{index} of i corresponds to the model with tuning parameters \code{object$all.alpha[i]} and \code{object$all.lambda[i]}.
}
}
\details{
The degrees of freedom for FLAM were derived in Section 4.1 of Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\value{
The degrees of freedom for the specified model.
}
\references{
Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\author{
Ashley Petersen
}

\examples{
#See ?'flam-package' for a full example of how to use this package

#generate data
#note: use larger 'n' for more reasonable results
set.seed(1)
data <- sim.data(n = 20, scenario = 1, zerof = 10, noise = 1)

#fit model for a range of tuning parameters
flam.out <- flam(x = data$x, y = data$y)
#or fit model and select tuning parameters using 2-fold cross-validation
#note: use larger 'n.fold' (e.g., 10) in practice
flamCV.out <- flamCV(x = data$x, y = data$y, n.fold = 2)

#calculate degrees of freedom for the model chosen using cross-validation
flamDOF(object = flamCV.out$flam.out, index = flamCV.out$index.cv)
#or for any fit from a 'flam' object
flamDOF(object = flam.out, index = 25)
flamDOF(object = flamCV.out$flam.out, index = 25)
#which corresponds to lambda and alpha of
flam.out$all.lambda[25]; flam.out$all.alpha[25]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
