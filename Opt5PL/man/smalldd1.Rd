
\name{smalldd1}
\alias{smalldd1}
\title{Sub-function of the function D_weight_2}
\usage{
smalldd1(T, x1, x2, xl, inv, order)
}
\arguments{
\item{T}{A numeric vector. Model parameter values}

\item{x1}{A numeric. ith dose level}

\item{x2}{A numeric. jth dose level}

\item{xl}{A numeric. last dose level}

\item{inv}{A numeric matrix. Inverse information matrix}

\item{order}{numeric. The number of model paraemters: 3, 4, 5 for the 3, 4, 5PL models respectively}
}
\description{
This computes the second derivative of the D-optimality criterion with respect to the model parameters for a specified model.
}
