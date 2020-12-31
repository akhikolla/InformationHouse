\name{dfphase1-package}
\alias{dfphase1-package}
\alias{dfphase1}
\docType{package}
\title{
  \packageTitle{dfphase1}
}
\description{
  \packageDescription{dfphase1}
}
\details{
  The main functions are:
  \itemize{
    \item \code{\link{shewhart}} and \code{\link{mshewhart}}:
    univariate and multivariate Shewhart-type control charts
    based either on the original observations or on a rank transformation.
    These functions are particularly useful for detecting isolated shifts
    in the mean and/or variance of subgrouped observations.
    Functions \code{\link{shewhart}} and \code{\link{mshewhart}} also allow
    the simultaneously use of two control charts originally
    designed to detect separately location and scale shifts.
    In particular, note that when more than one critical values are needed, the
    false alarm probability is \dQuote{balanced} between the
    separate control charts as discussed by Capizzi (2015).
    \item \code{\link{changepoint}} and \code{\link{mchangepoint}}:
    univariate or multivariate control charts useful for detecting
    sustained (and other patterned) mean and/or variance shifts.
    The control statistic is based on a generalized likelihood
    ratio test computed under a Gaussian assumption. However, the
    control limits are computed by permutation. An optional 
    preliminary rank transformation can be used to improve the performance in the case of
    nonnormal data.
    \item
    \code{\link{rsp}} and \code{\link{mphase1}}: the univariate and
    multivariate methods introduced by
    Capizzi and Masarotto (2013) and (2016) to detect multiple isolated
    or step shifts in individual or subgrouped data.
  }
  The use of distribution-free control limits is emphasized.
  However, the package also includes some functions 
  for computing normal-based control limits.
  As noted in the individual help pages, these limits can also be
  suitable  for some non-normal distributions
  (e.g., applying a multivariate rank.-transformation, normal-based
  control limits mantain the desired false alarm probability 
  in the class of the multivariate elliptical distributions).
  Nevertheless, their use is not generally recommended.

  The data should be organized as follows:
  \itemize{
    \item \emph{Univariate control charts:}
    an nxm matrix, where
    n and m are the size of each subgroup and
    the number of subgroups, respectively.
    A vector of length m is accepted  in the case of individual data, i.e.,
    when n=1.
    \item \emph{Multivariate control charts:}
    a pxnxm array, where
    p denotes the number of monitored variables.
    A p x m matrix is accepted in the case of individual
    data.
  }
  Functions \code{\link{phase1Plot}} and \code{\link{mphase1Plot}} can
  be used for plotting the data.
}
\author{
\packageAuthor{dfphase1} (maintainer: \packageMaintainer{dfphase1}).
}
\references{
  G. Capizzi (2015) \dQuote{Recent advances in process monitoring: Nonparametric and variable-selection
    methods for Phase I and Phase II (with discussion)}. \emph{Quality
    Engineering},  \bold{27}, pp. 44--80.
  
  G. Capizzi and G. Masarotto (2013),
  \dQuote{Phase I Distribution-Free Analysis of Univariate Data}.
  \emph{Journal of Quality Technology}, \bold{45}, pp. 273-284.
  
  G. Capizzi and G. Masarotto (2016),
  \dQuote{Phase I Distribution-Free Analysis with the R
    Package \code{dfphase1}}. To appear on
  \emph{Frontiers in Statistical Quality Control 12}, eds. S. Knoth and
  W. Schmid, Springer.

  G. Capizzi and G. Masarotto (2016), Phase I Distribution-Free
  Analysis of Multivariate Data, \emph{Technometrics},
  to appear, \url{http://dx.doi.org/10.1080/00401706.2016.1272494}.
  
}
\seealso{
  \code{\link{shewhart}},
  \code{\link{shewhart.normal.limits}},
  \code{\link{mshewhart}},
  \code{\link{mshewhart.normal.limits}},
  \code{\link{changepoint}},
  \code{\link{changepoint.normal.limits}},
  \code{\link{mchangepoint}},
  \code{\link{mchangepoint.normal.limits}},
  \code{\link{rsp}},
  \code{\link{phase1Plot}},
  \code{\link{mphase1Plot}}.
}
\keyword{ package }


