\name{clogitboost}
\alias{clogitboost}
\title{Boosting conditional logit model
}

\description{Fit a boosting conditional logit model using componentwise smoothing spline.
}

\usage{
clogitboost(y, x, strata, iter, rho)
}

\arguments{
  \item{y}{vector of binary outcomes.
}
  \item{x}{matrix or data frame with each column being a covariate.}
 \item{strata}{vector of group membership, i.e., items in the same group have the same value.}
 \item{iter}{number of iterations.}
 \item{rho}{learning rate parameter in the boosting algorithm.}
}

\value{
The function \code{clogitboost} returns the following list of values:
\item{call}{original function call.}
\item{func}{list of fitted spline functions.}
\item{index}{list of indices indicating which covariate is used as input for the smoothing spline.}
\item{theta}{list of fitted coefficients in the conditional logit models.}
\item{loglike}{sequence of fitted values of log-likelihood.}
\item{infscore}{relative influence score for each covariate.}
\item{rho}{learning rate parameter, which typically takes a value of 0.05 or 0.1.}
\item{xmax}{maximal element of each covariate.}
\item{xmin}{minimal element of each covariate.}

}


\author{
Haolun Shi \email{shl2003@connect.hku.hk}

Guosheng Yin \email{gyin@hku.hk}
}

\seealso{
\code{\link{plot.clogitboost}}

\code{\link{predict.clogitboost}}
}
\examples{
data(travel)
train <- 1:504
y <- travel$MODE[train]
x <- travel[train, 3:6]
strata <- travel$Group[train]
fit <- clogitboost(y = y, x = x, strata = strata, iter = 10, rho = 0.05)
}

