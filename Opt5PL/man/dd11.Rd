
\name{dd11}
\alias{dd11}
\title{Computing each element of the function DD_weight_2}
\usage{
dd11(T, x1, x2, xl, inv, inv1, order)
}
\arguments{
\item{T}{A numeric vector. Model parameter values}

\item{x1}{A numeric. ith dose level for a given design}

\item{x2}{A numeric. jth dose level for a given design}

\item{xl}{A numeric. last dose level for a given design}

\item{inv}{A numeric matrix. Inverse of the information matrix for the 5PL model}

\item{inv1}{A numeric matrix. Inverse of the information matrix for the 4PL model}

\item{order}{numeric. The number of model paraemters}
}
\description{
This function computes each element of the matrix function DD_weight_2 which is the second derivative of the Ds-optimality criterion with respect to the model parameters.
}
