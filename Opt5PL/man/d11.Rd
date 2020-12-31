
\name{d11}
\alias{d11}
\title{Computing each element of the function DD_weight_1}
\usage{
d11(T, x, xl, inv, inv1, order)
}
\arguments{
\item{T}{A numeric vector. Model parameter values}

\item{x}{A numeric. ith dose level for a given design}

\item{xl}{A numeric. last dose level for a given design}

\item{inv}{A numeric matrix. Inverse of the information matrix for the 5PL model}

\item{inv1}{A numeric matrix. Inverse of the information matrix for the 4PL model}

\item{order}{numeric. The number of model paraemters}
}
\description{
This function computes each element of the vector function DD_weight_1 which is the first derivative of the Ds-optimality criterion with respect to the model parameters.
}
