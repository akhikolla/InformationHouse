\name{addRowWiseC}
\alias{addRowWiseC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
intern function
}
\description{
 Adds the Vector DataPoint to every row of the matrix WeightVectors}
\usage{
addRowWiseC(WeightVectors,DataPoint)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{WeightVectors}{ WeightVectors. n weights with m components each}
    \item{DataPoint}{Vector with m components}
}

\value{
WeightVectors[1:m,1:n]
}