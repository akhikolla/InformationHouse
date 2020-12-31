\name{HWQqplot}
\alias{HWQqplot}
\alias{HWQqplot}
\title{
A Q-Q plot for Hardy-Weinberg equilibrium
}
\description{
\code{HWQqplot} creates a Q-Q plot for the p-values obtained in an Exact
test for Hardy-Weinberg equilibrium. Empirical p-values are plotted
against multiple simulated quantiles of the theoretical p-value distribution.
}
\usage{
HWQqplot(X, nsim = 100, fit = "curve", logplot = FALSE,
main = "Q-Q plot for HWE", mm = NULL, pvaluetype = "selome", ...)
}
\arguments{
  \item{X}{Data matrix with genotype counts, one row for each sample, 3 columns}
  \item{nsim}{Number of samples drawn from the null distribution (100 by
  default)}
  \item{fit}{If \code{fit} is set to "line" straight lines will be
    fitted to the simulated samples, if set to "curve", ascending curves
    will be shown.}
  \item{logplot}{If \code{logplot} is set to true, then the log10 of the
  p-values will be used in the plot. If not, untransformed p-values will
  be used.}
  \item{main}{Title for the plot}
  \item{mm}{Maximal value for x and y axis in the plot}
  \item{pvaluetype}{Type of p-value to be used in an exact test. Can be
    "selome" (default), "midp" or "dost".}
  \item{\dots}{Any additional arguments for the \code{plot} instruction}
}
\details{
\code{HWQqplot} constructs a Q-Q plot of the p-values of an exact test
for Hardy-Weinberg equilibrium. Under the null, this p-value is not
uniform. HWQqplot samples from the theoretical null distribution, taking
into account that markers may vary in allele frequency and in sample
size (due to missing values). For each simulated sample a grey curve or
line is shown. A green reference line with intercept 0 and slope 1 is
also shown in the plot.
}
\value{NULL}
\references{
Rohlfs, R.V. and Weir, B.S. (2008) Distributions of Hardy-Weinberg
equilibrium test statistics. Genetics 180, pp. 1609-1616.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWTernaryPlot}} \code{\link{HWExact}} \code{\link{qqplot}}}
\examples{
\dontrun{
set.seed(1234)
n <- 200 # sample size
m <- 100 # number of markers
X <- HWData(n,m)$Xt
HWQqplot(X,logplot=TRUE,pvaluetype="selome",main="Q-Q Plot for HWE")
}
}
\keyword{aplot}

