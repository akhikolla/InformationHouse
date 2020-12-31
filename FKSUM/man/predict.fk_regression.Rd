\name{predict.fk_regression}
\alias{predict.fk_regression}
\title{Predict method for class fk_regression}
\description{
  Predict method for kernel regression. Evaluates the fitted regression function at a set of evaluation/test points
}
\usage{
\method{predict}{fk_regression}(object, xtest = NULL, ...)
}
\arguments{
  \item{object}{an object of class fk_regression, output from the function of the same name.}
  \item{xtest}{(optional) vector of evaluation/test points. If omitted then the fitted values at the original sample are returned.}
  \item{...}{(optional) further arguments passed to or from other methods.}
}
\value{vector of fitted/predicted values of the regression function.}
\keyword{file}
