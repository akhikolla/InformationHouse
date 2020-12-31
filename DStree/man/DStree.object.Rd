\name{DStree.object}
\alias{DStree.object}
\title{
  DStree Object 
}
\description{
  This class of objects is returned by the \code{\link{DStree}} function
  to represent a fitted decision tree for right-censored survival data.}
\section{Structure}{
  The above components must be included in a legitimate \code{DStree} object.
}
\value{

\item{frame}{
  data frame with one row for each node in the tree.
  The \code{row.names} of \code{frame} contain the (unique) node numbers that
  follow a binary ordering indexed by node depth.
  Columns of \code{frame} include
  \code{var}, a factor giving the names of the variables used in the
  split at each node (leaf nodes are denoted by the level \code{"<leaf>"}), 
  \code{n}, the number of observations reaching the node,
  \code{wt}, the sum of case weights for observations reaching the node,
  \code{dev}, the deviance of the node,
  \code{yval}, the estimated value of the median survival time at the node,
  and \code{splits}, a two-column matrix of left and right split labels
  for each node. Also included in the frame are \code{complexity}, the
  complexity parameter at which this split will collapse, \code{ncompete},
  the number of competitor splits recorded, and \code{nsurrogate}, the
  number of surrogate splits recorded. \code{yval2} contains in its first columns the fitted hazard probabilities and in last the fitted survival probabilities (see \code{\link{surv}}).
}
\item{where}{
  an integer vector of the same length as the number of observations in the
  root node, containing the row number of \code{frame} corresponding to
  the leaf node that each observation falls into. 
}
\item{call}{
  an image of the call that produced the object, but with the arguments 
  all named and with the actual formula included as the formula argument. 
  To re-evaluate the call, say \code{update(tree)}.
}
\item{terms}{
  an object of class \code{c("terms", "formula")} (see
  \code{\link{terms.object}}) summarizing the formula.  Used by various
  methods, but typically not of direct relevance to users.
}
\item{splits}{
  a numeric matrix describing the splits: only present if there are any.
  The row label is the name of
  the split variable, and columns are \code{count}, the number of
  observations (which are not missing and are of positive weight) sent
  left or right by the split (for competitor splits this is the number
  that would have been sent left or right had this split been used, for
  surrogate splits it is the number missing the primary split variable
  which were decided using this surrogate), \code{ncat}, the number of
  categories or levels for the variable (\code{+/-1} for a continuous
  variable), \code{improve}, which is the improvement in deviance given
  by this split, or, for surrogates, the concordance of the surrogate
  with the primary, and \code{index}, the numeric split point.  The last
  column \code{adj} gives the adjusted concordance for surrogate splits.
  For a factor, the \code{index} column contains the row number of the
  csplit matrix.  For a continuous variable, the sign of \code{ncat}
  determines whether the subset \code{x < cutpoint} or \code{x >
  cutpoint} is sent to the left.
}
\item{csplit}{
  an integer matrix.  (Only present only if at least one of the split
  variables is a factor or ordered factor.)  There is a row for
  each such split, and the number of columns is the largest number of
  levels in the factors.  Which row is given by the \code{index} column
  of the \code{splits} matrix.  The columns record \code{1} if that
  level of the factor goes to the left, \code{3} if it goes to the
  right, and \code{2} if that level is not present at this node
  of the tree (or not defined for the factor).
}
\item{method}{
  character string: the method used to grow the tree. Since \code{DStree} is based on a user-defined split function of \code{rpart} the method is always 
  \code{"user"}.
}
\item{cptable}{
  a matrix of information on the optimal prunings based on a
  complexity parameter.
}
\item{variable.importance}{
  a named numeric vector giving the importance of each variable.  (Only
  present if there are any splits.)  When printed by
  \code{\link{summary.DStree}} these are rescaled to add to 100.
}
\item{numresp}{
  integer number of responses; the number of levels for a factor response.
}
\item{parms, control}{
  a record of the arguments supplied, which defaults filled in.
}
\item{functions}{
  the \code{summary} and split function used.
}
\item{ordered}{
  a named logical vector recording for each variable if it was an
  ordered factor.
}
\item{na.action}{
  (where relevant) information returned by \code{\link{model.frame}} on
  the special handling of \code{NA}s derived from the \code{na.action}
  argument.
}
\item{names}{a string vector of size two, which denotes the column names of the observed time points and of the \code{status} variable.

}
\item{wt}{
a numeric vector of equal length as the number of rows of the dataset, which denotes the optional case weigths, defined in the \code{weights} argument.
}


Optional components include the model frame (\code{model}), the matrix
of predictors (\code{x}) and the response variable (\code{y}) used to
construct the \code{DStree} object.
}
\seealso{
  \code{\link{DStree}}.
}
\keyword{tree}
\keyword{methods}
