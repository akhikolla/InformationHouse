\name{evalStepFun}
\alias{evalStepFun}

\title{Evaluate step function}
\description{Transform the return value by \code{\link{smuce}}, \code{\link{fdrseg}}, or \code{\link{dfdrseg}} into a numeric vector.}

\usage{
evalStepFun(stepF)
}

\arguments{
\item{stepF}{a list returned by \code{\link{smuce}}, \code{\link{fdrseg}}, or \code{\link{dfdrseg}} with components 
\preformatted{
value function values on each segment of the estimator
left  indices of leftmost points within each segment of the estimator
n     number of samples
}
} 
}

\value{
A numeric vector gives function values of \code{stepF} at sampling locations.  
}
 

\seealso{\code{\link{smuce}}, \code{\link{fdrseg}}, \code{\link{dfdrseg}}}

\examples{
# simulate data
set.seed(2)
u0 <- c(rep(1, 5), rep(5, 5))
Y  <- rnorm(10, u0)

# compute the SMUCE estimate
uh <- smuce(Y)

# print results
#   step function returned by smuce
print(uh)
#   vector returned by evalStepFun
print(evalStepFun(uh))
}

\keyword{nonparametric}


