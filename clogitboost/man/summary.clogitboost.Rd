\name{summary.clogitboost}
\alias{summary.clogitboost}

\title{
Summary after fitting a boosting conditional logit model
}

\description{
  \code{summary} methods for the \code{clogitboost} objects.
}
  
\usage{
\method{summary}{clogitboost}(object, \dots)
}

\arguments{
    \item{object}{output object from the \code{\link{clogitboost}} function.
        }
    
    \item{\dots}{not currently used.}
}



\value{
The function \code{clogitboost()} returns the following list of values:
\item{call}{original function call.}
\item{infscore}{relative influence score for each covariate.}
\item{loglike}{sequence of the fitted values of log-likelihood.}


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
summary(fit)
}

