\name{DStree.control}
\alias{DStree.control}
\title{Control arguments for DStree algorithm}
\usage{
DStree.control(minsplit = 20L, minbucket = round(minsplit/3), cp = 0.005,
  maxcompete = 4L, maxdepth = 30, maxsurrogate = 0)
}
\arguments{
  \item{minsplit}{the minimum number of observations that
  must exist in a node in order for a split to be
  attempted.}

  \item{minbucket}{the minimum number of observations in
  any terminal node. If only one of \code{minbucket}
  or \code{minsplit} is specified, the code either sets
  \code{minsplit} to \code{minbucket}*3 or \code{minbucket}
  to \code{minsplit}/3, as appropriate.}

  \item{cp}{complexity parameter. Any split that does not
  decrease the overall lack of fit by a factor of cp is not
  attempted. This means that the overall deviance must
  decrease by cp at each step. The main role of this
  parameter is to save computing time by pruning off splits
  that are not worthwhile by definition.  Essentially, the user
  informs the program that any split which does not improve
  the fit by cp will likely be pruned off.}

  \item{maxcompete}{the number of competitor splits
  retained in the output. It is useful to know not just
  which split was chosen, but which variable came in
  second, third, etc.}

  \item{maxdepth}{Set the maximum depth of any node of the
  final tree, with the root node counted as depth 0. Values
  greater than 30 \code{DStree} will give nonsense results
  on 32-bit machines.}

  \item{maxsurrogate}{the number of surrogate splits
  retained in the output. If this is set to zero the
  compute time will be reduced, since approximately half of
  the computational time (other than setup) is used in the
  search for surrogate splits.}
}

\description{
Various paramters that control the specification of the \code{DStree}
algorithm.
}

