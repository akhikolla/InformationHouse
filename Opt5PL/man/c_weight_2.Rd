
\name{c_weight_2}
\alias{c_weight_2}
\title{The second derivative of the c-optimality criterion with respect to the model parameters}
\usage{
c_weight_2(W, T, X, inv, p, order)
}
\arguments{
\item{W}{A numeric vector.  The first K-1 weights for a given design}

\item{T}{A numeric vector. Model parameter values}

\item{X}{A numeric vector. K dose levels for a given design}

\item{inv}{A numeric matrix. Inverse of the information matrix.}

\item{p}{A numeric. Define EDp}

\item{order}{numeric. The number of model paraemters}
}
\description{
The second derivative of the c-optimality criterion with respect to the model parameters. This is a sub-function of c_weight to run one iteration of Newton-Raphson method.
}
