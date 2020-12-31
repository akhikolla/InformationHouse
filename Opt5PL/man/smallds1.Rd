
\name{smallds1}
\alias{smallds1}
\title{Sensitivity function of D-optimality criterion}
\usage{
smallds1(T, x, inv, order)
}
\arguments{
\item{T}{A numeric vector. Model parameter values}

\item{x}{A numeric. A given dose level}

\item{inv}{A numeric matrix. Inverse information matrix}

\item{order}{numeric. The number of model paraemters: 3, 4, 5 for the 3, 4, 5PL models respectively}
}
\description{
This evaluates the sensitivity function of the D-optimality criterion at the given value of x
 for a specified model.}
