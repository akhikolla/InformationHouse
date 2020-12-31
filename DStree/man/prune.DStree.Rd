\name{prune.DStree}
\alias{prune.DStree}
\title{Prune a fitted survival tree}
\usage{
\method{prune}{DStree}(tree, data, gamma = 2, which, ...)
}
\arguments{
  \item{tree}{fitted model of class \code{"DStree"}.
  This object is assumed to be the result of some function that
  produces an object with the same named components as that
  returned by the \code{DStree} function.}

  \item{data}{optional data frame that is used to evaluate
  the fit of the tree. The
  predictors referred to in \code{formula(tree)}
  as well as the 'status' variable must be present by name
  in \code{data}. If missing, the data frame that was used to
  grow the tree is used.}

  \item{gamma}{optional positive integer value that determines
  the user defined information criterion. Setting \code{gamma} = 2
  corresponds to the AIC criterium.}

  \item{which}{An optional string that determines which performance criteria
  should be computed from the subtrees. One of \code{"CRIT"}, \code{"DEV"}
  or \code{"BRIER"}. If \code{which}
  is missing the \code{"BRIER"} score is chosen.}

  \item{...}{further arguments passed to or from other
  methods.}
}
\value{
\code{prune} returns one \code{DStree} object and four vectors of length equal to the
number of subtrees: \itemize{
\item \code{nsplit} number of splits for every subtree
\item \code{CRIT} value of the user defined information criterion (underlying
formula: \emph{CRIT = deviance + gamma * |terminal leaves| * |time
periods|}). \item \code{DEV} deviance \item \code{BRIER}
Integrated Brier Score, see Hothorn et al. (2004) \item
\code{pruned.fit} optimal subtree regarding the choosen
criterium specified in \code{which}}
}
\description{
The \code{prune} function evaluates and prunes a survival tree that has
been fitted by \code{DStree}. Different criteria can be used for evaluation
(e.g. Brier Score or AIC). The result of the \code{prune} function
is the optimal subtree (of class \code{DStree}) with regard to the chosen
criterium, as well as various performance measures that were obtained from the subtrees
during analysis. The returned performance measures are the Brier Score, the deviance,
and an information criterion defined by \code{gamma}.
}
\details{
The subtrees are the cost-minimzing subtrees in terms of
deviance for given complexity parameters of the fitted
tree. See Therneau et al (2013) p.12-13.
}
\examples{
data(cost)
## Discretize observed days to years
d.cost <- dis.cost(cost)

##Build tree
tree <- DStree(time~prevStroke+age+sex+alcohol+smoke,status="status",data=d.cost)

# Determine subtree with minimum AIC
prunedtree <- prune(tree,d.cost,which="CRIT")
prunedtree$prunedfit

# Visualize AIC/Deviance of subtrees
plot(prunedtree$nsplit,prunedtree$CRIT)
plot(prunedtree$nsplit,prunedtree$DEV)
}
\references{
Hothorn T., Lausen B., Benner A. and Radespiel-Troeger M.
(2004), Bagging Survival Trees. \emph{Statistics in
medicine 23 (1), 77-91}.

Therneau T. and Atkinson E., An
introduction to recursive partitioning using the RPART
routines, \emph{Technical Report 61, Section of
Biostatistics, Mayo Clinic, Rochester.}
}

