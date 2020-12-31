
\name{DD_weight_2}
\alias{DD_weight_2}
\title{The second derivative of the Ds-optimality criterion with respect to the model parameters}
\usage{
DD_weight_2(W, T, X, inv, inv1, order)
}
\arguments{
\item{W}{A numeric vector. The first K-1 weights for a given design}

\item{T}{A numeric vector. Model parameter values for the 5PL model}

\item{X}{A numeric vector.  K dose levels for a given design}

\item{inv}{A numeric matrix. Inverse information matrix for the 5PL model}

\item{inv1}{A numeric matrix. Inverse information matrix for the 4PL model}

\item{order}{numeric. The number of model paraemters}
}
\description{
The second derivative of the Ds-optimality criterion with respect to the model parameters. This is a sub-function of D_weight to run one iteration of Newton-Raphson method.
}
