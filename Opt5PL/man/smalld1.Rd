
\name{smalld1}
\alias{smalld1}
\title{Sub-function of the function D_weight_1}
\usage{
smalld1(T, x, xl, inv, order)
}
\arguments{
\item{T}{A numeric vector. Model parameter values}

\item{x}{A numeric. ith dose level}

\item{xl}{A numeric. last dose level}

\item{inv}{A numeric matrix. Inverse infromation matrix}

\item{order}{numeric. The number of model paraemters: 3, 4, 5 for the 3, 4, 5PL models respectively}
}
\description{
This computes the first derivative of the D-optimality criterion with respect to the model parameters for a specified model.
}
