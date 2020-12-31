\name{predict.clogitboost}
\alias{predict.clogitboost}

\title{Predicting after fitting a boosting conditional logit model}

\description{
  \code{predict} methods for the \code{clogitboost} objects, which produce marginal predictions of the covariate effects.
}
  
\usage{
\method{predict}{clogitboost}(object, x, strata, \dots)
}

\arguments{
    \item{object}{output object from the \code{\link{clogitboost}} function.
        }
     \item{x}{new matrix or data frame with each column being a covariate.}
 \item{strata}{new vector of group memberships, i.e., items in the same group have the same value.}
    \item{\dots}{not currently used.}
}


\value{
The method \code{predict} returns the following list of values:
\item{prob}{probability of the outcome equal to 1.}
\item{utility}{predicted utility.}
\item{prediction}{0-1 prediction of the outcome variable.}


}

\author{
Haolun Shi \email{shl2003@connect.hku.hk}

Guosheng Yin \email{gyin@hku.hk}
}

\seealso{
\code{\link{clogitboost}}
}
\examples{
data(travel)
train <- 1:504
y <- travel$MODE[train]
x <- travel[train, 3:6]
strata <- travel$Group[train]
fit <- clogitboost(y = y, x = x, strata = strata, iter = 10, rho = 0.05)
predict(fit, x = travel[-train, 3:6], strata = travel$Group[-train])
}

