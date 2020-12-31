\name{IsingPL}
\alias{IsingPL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudo-likelihood
}
\description{
Computes the pseudolikelihood of a dataset given an Ising Model.
}
\usage{
IsingPL(x, graph, thresholds, beta, responses = c(0L, 1L))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A binary dataset
}
  \item{graph}{
Square matrix indicating the weights of the network. Must be symmetrical with 0 as diagonal.
}
  \item{thresholds}{
Vector indicating the thresholds, also known as the external field.
}
  \item{beta}{
Scalar indicating the inverse temperature.
}
  \item{responses}{
Response options. Typically set to \code{c(-1L, 1L)} or \code{c(0L, 1L)} (default). Must be integers!
}
}

\value{
The pseudolikelihood
}
\author{
Sacha Epskamp (mail@sachaepskamp.com)
}
