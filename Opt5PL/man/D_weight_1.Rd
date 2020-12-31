
\name{D_weight_1}
\alias{D_weight_1}
\title{The first derivative of the D-optimality criterion w.r.t the model parameters}
\usage{
D_weight_1(q, W, T1, T2, T3, X, inv1, inv2, inv3)
}
\arguments{
\item{q}{A numeric vector. User select weights defined in the function RDOPT}

\item{W}{A numeric vector. The first K-1 weights for a given design}

\item{T1}{A numeric vector. Model parameter values for the 3PL model}

\item{T2}{A numeric vector. Model parameter values for the 4PL model}

\item{T3}{A numeric vector. Model parameter values for the 5PL model}

\item{X}{A numeric vector. K dose levels for a given design}

\item{inv1}{A numeric matrix. Inverse information matrix for the 3PL model}

\item{inv2}{A numeric matrix. Inverse information matrix for the 4PL model}

\item{inv3}{A numeric matrix. Inverse information matrix for the 5PL model}
}
\description{
The first derivative of the D-optimality criterion with respect to the model parameters. This is a sub-function of D_weight to run one iteration of Newton-Raphson method.
}
