\name{make.tree}
\alias{make.tree}

\title{
Growing a binary tree for change-point analysis
}

\description{
Grows a binary tree of a given height for input panel data,
with possible changes in the means of component series,
using Double CUSUM statistics, where each node contains information
about estimated change-points
}

\usage{
make.tree(y, dw, rule = NULL)
}

\arguments{
  \item{y}{
  a matrix containing each component series as its rows
  }
  \item{dw}{
  see \code{dw} in \code{\link{factor.seg.alg}}
  }
  \item{rule}{
  see \code{rule} in \code{\link{factor.seg.alg}}
  }
}

\value{
  \item{tree}{
  a list representing a binary tree grown as an outcome from
  the Double CUSUM Binary Segmentation algorithm;
  each node is a \code{6}-by-\code{1} matrix containing
  the index of the node, the start and end of the segment and
  an estimated change-point candidate, and the test statistic
  }
  \item{mat}{
  a matrix stacking the nodes of \code{tree} as its columns
  }
}

\author{
Haeran Cho
}

\references{
M. Barigozzi, H. Cho and P. Fryzlewicz (2016) Simultaneous multiple change-point and factor analysis for high-dimensional time series, Preprint.
}

\keyword{internal}
