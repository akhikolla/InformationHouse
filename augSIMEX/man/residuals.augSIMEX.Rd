\name{residuals.augSIMEX}
\alias{residuals.augSIMEX}

\title{
Residuals of the Fits Made by augSIMEX
}
\description{
This function outputs the residuals from a fitted augSIMEX object.
}
\usage{
## S3 method for class 'augSIMEX'

\method{residuals}{augSIMEX}(object, type = c("deviance", "pearson", "working",
   "response", "partial"), ...) 
}

\arguments{
  \item{object}{
the ``augSIMEX" object gotten from \code{\link[augSIMEX]{augSIMEX}} function.
}
  \item{type}{
a character specifying what type of residuals to be returned. The option includes the deviance, Pearson, working, response and partial.
}
  \item{\dots}{
other arguments that are passed into the function.
}
}

\author{
Qihuang Zhang and Grace Y. Yi.
}

\seealso{
\code{\link[stats]{residuals}},\code{\link[stats]{residuals.glm}}
}

\examples{
data(ToyUni)
example <- augSIMEX(mainformula = Y ~ Xstar + Zstar + W, family = binomial(link = logit),
  mismodel = pi|qi ~ W, 
  meformula = Xstar ~ X + Z + W,
  data = ToyUni$Main,validationdata = ToyUni$Validation, subset = NULL,
  err.var = "Xstar", mis.var = "Zstar", err.true = "X", mis.true = "Z", 
  err.mat = NULL,
  lambda = NULL, M = 5, B = 2, nBoot = 2, extrapolation="quadratic")
residuals(example)
}

\keyword{function}
