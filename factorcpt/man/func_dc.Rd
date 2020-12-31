\name{func_dc}
\alias{func_dc}

\title{
Double CUSUM statistics
}

\description{
Aggregates high-dimensional CUSUM statistic series via Double CUSUM approach.
}

\usage{
func_dc(z)
}

\arguments{
  \item{z}{
  input data matrix
  }
}

\value{
\item{cs}{a matrix of CUSUM series}
\item{acs}{a matrix of absolute values of CUSUM series}
\item{res}{pointwise maximum of Double CUSUM statistics}
\item{mat}{a matrix of Double CUSUM statistics}
}

\references{
H. Cho (2016) Change-point detection in panel data via double CUSUM statistic. Electronic Journal of Statistics. 10: 2000-2038.}

\author{
Haeran Cho
}

\examples{
e <- matrix(rnorm(100*100), nrow=100)
f <- matrix(0, nrow=100, ncol=100)
f[, 51:100] <- 1 # t = 50 represents the change-point in the means
x <- f + e
fd <- func_dc(x)
plot(fd$res, type='l'); abline(v=50, col=2) # pointwise maximum of DC statistics
image(t(fd$mat)) # heatmap of DC statistics with x-axis representing the time
}
