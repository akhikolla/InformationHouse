\name{vs.test}
\alias{vs.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vasicek-Song goodness-of-fit test for various distributions}
\description{Performs Vasicek-Song goodness-of-fit test to the specified distribution family.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vs.test(x, densfun, param = NULL, 
        simulate.p.value = NULL, B = 5000,
        delta = NULL, extend = FALSE, relax = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{(\code{numeric, vector}) the numeric sample.}
  \item{densfun}{A character string specifying the fitted distribution. Possible values are \code{"dunif"}, \code{"dnorm"}, \code{"dlnorm"}, \code{"dexp"}, \code{"dgamma"}, \code{"dweibull"}, \code{"dpareto"}, \code{"df"}, \code{"dlaplace"} and \code{"dbeta"}.}
  \item{param}{(\code{numeric, vector}) specifies the parameter(s) of the fitted distribution. If \code{NULL} (default), a GOF test to the parametric family of distributions specified by densfun is performed.}
  \item{simulate.p.value}{(\code{logical, single value}) if \code{TRUE}, the p-value of the sample is estimated by means of Monte Carlo methods. If \code{NULL} (the default), the p-value is simulated if the sample size is smaller than 80; otherwise, an asymptotic p-value is computed.}
  \item{B}{(\code{numeric, single value}) a numeric value specifying the number of simulations to perform in Monte-Carlo estimation of the p-value.}
  \item{delta}{(\code{numeric, single value}) a numeric value smaller than \eqn{1/3} specifying the upper bound \eqn{n^{1/3}-\delta} for window size, where \eqn{n} is the sample size. The default depends on \code{densfun}; see Vignettes for details.}
  \item{extend}{(\code{logical, single value}). If \code{FALSE} (the default), the bound  for the window is \eqn{n^{1/3}-\delta}; if \code{TRUE}, the bound is \eqn{n/2}.}
  \item{relax}{(\code{logical, single value}) avoids the constraint \eqn{V_{mn} \leq -\frac{1}{n} \sum_{i=1}^n \log p_0(X_i, \widehat{\theta}_n)} when computing the optimal window; see details. Default is \code{FALSE}.}
}

\details{The test statistic is 
\deqn{I_{mn}=-V_{mn}-\frac{1}{n}\sum_{i=1}^{n}\log p_{0}(X_{i},\theta),}
where \eqn{V_{mn}} is the Vasicek estimator of Shannon entropy computed from the numeric sample \code{x} with window size \eqn{m} and \eqn{p_{0}(x,\theta)} is the density function of the specified distribution \code{densfun} to be tested, with \eqn{\theta} the parameter of the null for a simple hypothesis or its maximum likelihood estimate for a composite null hypothesis (\code{param=NULL}); See Song (2002), Girardin and Lequesne (2017) and Lequesne and Regnault (2018). 

An optimal window size \eqn{m} is automatically computed; see Song (2002).

An  exact p-value is computed if the sample size is less than 100. Otherwise, asymptotic distribution is used whose approximation may be inaccurate for small samples; see Lequesne and Regnault (2018).
}
\value{ A list with class "htest" containing the following components:
  \item{observed}{The sample under study.}
  \item{data.name}{The name (as an R object) of the sample.}
  \item{null.value}{A character string specifying the name of the fitted distribution.}
  \item{method}{The character string \code{"Vasicek GOF test to"} followed by the name of the fitted distribution.}
  \item{statistic}{Vasicek test statistic; see Details below.}
  \item{parameter}{The optimal window for Vasicek test statistic}
  \item{estimate}{Parameter(s) of the fitted distribution. If \code{param} is \code{NULL}, parameters are estimated. If \code{param} is suitably filled out by the user, it is returned.}
  \item{p.value}{The p-value of the test.}
}
\references{
Vasicek, O., A test for normality based on sample entropy, \emph{Journal of the Royal Statistical Society,} \bold{38(1)}, 54-59 (1976).

Song, K. S., Goodness-of-fit tests based on Kullback-Leibler discrimination information, \emph{Information Theory, IEEE Transactions on,} \bold{48(5)}, 1103-1117 (2002).

Girardin, V., Lequesne, J. Entropy-based goodness-of-fit tests - a unifying framework. Application to DNA replication. \emph{Communications in Statistics: Theory and Methods} (2017). https://doi.org/10.1080/03610926.2017.1401084

Lequesne, J., Regnault, P. vsgoftest: An R Package for Goodness-of-Fit Testing Based on Kullback-Leibler Divergence. \emph{Journal of Statistical Software}, \bold{96} (2020).  doi:10.18637/jss.v096.c01 

}
\author{
J. Lequesne \email{justine.lequesne@unicaen.fr}
}

\seealso{\code{\link{entropy.estimate}} which computes the Vasicek estimator of Shannon entropy.}

\examples{
set.seed(1)
samp <- rnorm(50,2,3)
vs.test(x = samp, densfun = 'dnorm', param = c(2,3), B = 500) #Simple null hypothesis
vs.test(x = samp, densfun='dnorm', B = 500) #Composite null hypothesis
## Using asymptotic distribution to compute the p-value
vs.test(x = samp, densfun='dnorm', simulate.p.value = FALSE) #Composite null hypothesis

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}