\name{plot.augSIMEX}
\alias{plot.augSIMEX}

\title{
Plot of Extrapolation
}
\description{
This function visualizes the extrapolation step. It plots the simulated value for each lambda value and the curve of extrapolation. The function provide the plots for both methods simultaneously, and guides the user to choose a proper extrapolation method. 
}
\usage{
## S3 method for class 'augSIMEX'

\method{plot}{augSIMEX}(x,...)
}

\arguments{
  \item{x}{
the ``augSIMEX" object gotten from \code{\link[augSIMEX]{augSIMEX}} function.
}
  \item{\dots}{
other arguments that are passed into the function.
}
}
\details{
The user may need to adjust the range of y axis for a proper display.
}
\author{
Qihuang Zhang and Grace Y. Yi.
}

\seealso{
\code{\link{plot}}
}

\examples{
data(ToyUni)
example<-augSIMEX(mainformula = Y ~ Xstar + Zstar + W, family = binomial(link = logit),
  mismodel = pi|qi ~ W, 
  meformula = Xstar ~ X + Z + W,
  data = ToyUni$Main,validationdata = ToyUni$Validation, subset = NULL,
  err.var = "Xstar", mis.var = "Zstar", err.true = "X", mis.true = "Z", 
  err.mat = NULL,
  lambda = NULL, M = 5, B = 2, nBoot = 2, extrapolation="quadratic")
         
plot(example, ylim = c(-1,0.4))
}

\keyword{function}
