\name{finity-package}
\alias{finity-package}
\alias{finity}
\docType{package}
\title{
\packageTitle{finity}
}
\description{
\packageDescription{finity}
}
\details{
The main function of the package, which performs the test for the finiteness of moments, is 

finite_moment_test().

The test assumes a sample drawn from an unknown distrubution \eqn{F}. Given this sample, the finiteness or infiniteness of any moment of  order \eqn{k} of distribution \eqn{F} can be ascertained. For this, the test follows a randomised testing procedure with artificial randomness. The absolute sample moment \eqn{\mu_k} of the desired order \eqn{k} of the sample is transformed into a test statistic,  which follows a \eqn{\chi^2} distribution with one degree of freedom exactly if the moment of the same order \eqn{k} of \eqn{F} is not finite. The null hypothesis in the test is that the moment is infinite; the alternative is that it is finite.

It should be noted that while the moment of order \eqn{k} of \eqn{F} may be infinite, the sample moment \eqn{\mu_k} is always finite because the sample is of finite size. The sample moment will, however, diverge with growing sample size if the moment of the same order \eqn{k} of the original distribution \eqn{F} is not finite.

The test works as follows: A standard normal distribution is rescaled with \eqn{\sqrt{\exp(\mu_k)}}, yielding a normal distribution with mean 0 and either finite or infinite variance, depending on whether the hypothesis holds. For every observation of the resulting distribution, it is then tested, if the observation is located within an interval \eqn{[-u, u]}. The resulting binary quantity \eqn{\zeta} (0, or 1, true or false) follows a Bernoulli distribution with mean 1/2 exactly if the \eqn{k}th moment of \eqn{F} is infinite. Otherwise the mean is not 1/2. Sampling a number of different intervals characterized by different bounds \eqn{u} drawn from a distribution with finite support, the test aggregates over quantities \eqn{\zeta} such that the resulting test statistic follows a \eqn{\chi^2} exactly if \eqn{E(\zeta)=1/2}, i.e., if the \eqn{k}th moment of \eqn{F} is infinite.

Trapani (2016) offers some insights into the performance of the test and the impact its parameters have. These parameters are optional arguments of the testing function finite_moment_test().
}
\author{
\packageAuthor{finity}

Maintainer: \packageMaintainer{finity}
}
\references{
Trapani, Lorenzo (2016), 'Testing for (in)finite moments', \emph{Journal of Econometrics} \bold{191}(1), 57 – 68.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\seealso{
https://github.com/x0range/finity
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# Generate sample
rvs <- stabledist::rstable(100000, 1.9, 0.5, 1, 0, pm = 0)
# Perform test
result <- finite_moment_test(rvs, 2)
# Print results
message(paste("Test statistic:", result[1], "p-value:", result[2], "\n\n"))

# More examples are included in https://github.com/x0range/finity/examples
}
