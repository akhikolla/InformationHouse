
\name{D1}
\alias{D1}
\title{Computing each element of the function c_weight_1}
\usage{
D1(T, x, xl, inv, p, order)
}
\arguments{
\item{T}{A numeric vector. Model parameter values}

\item{x}{A numeric. ith dose level for a given design}

\item{xl}{A numeric. last dose level for a given design}

\item{inv}{A numeric matrix. Inverse of the information matrix}

\item{p}{A numeric. Define EDp}

\item{order}{numeric. The number of model paraemters}
}
\description{
This function computes each element of the vector function c_weight_1 which is the first derivative of the c-optimality criterion with respect to the model parameters.
}
