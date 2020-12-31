
\name{DD1}
\alias{DD1}
\title{Computing each element of the function c_weight_2}
\usage{
DD1(T, x1, x2, xl, inv, p, order)
}
\arguments{
\item{T}{A numeric vector. Model parameter values}

\item{x1}{A numeric. ith dose level for a given design}

\item{x2}{A numeric. jth dose level for a given design}

\item{xl}{A numeric. last dose level for a given design}

\item{inv}{A numeric matrix. Inverse of the information matrix}

\item{p}{A numeric. Define EDp}

\item{order}{numeric. The number of model paraemters}
}
\description{
This function computes each element of the matrix function c_weight_2 which is the second derivative of the c-optimality criterion with respect to the model parameters.
}
