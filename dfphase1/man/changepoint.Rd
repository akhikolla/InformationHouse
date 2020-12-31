\name{changepoint}
\alias{changepoint}
\alias{mchangepoint}
\alias{changepoint.normal.limits}
\alias{mchangepoint.normal.limits}
\encoding{utf8}
\title{
  Detection of a sustained change-point in univariate and multivariate data
}
\description{
  \code{changepoint} (univariate data) and \code{mchangepoint}
  (multivariate data) test for the presence of a
  sustained location and/or dispersion shift. Both functions can be applied
  to individual and subgrouped observations.

  \code{changepoint.normal.limits} and
  \code{mchangepoint.normal.limits} precompute
  the corresponding control limits when the in-control distribution is
  normal.
}
\usage{
changepoint(x, subset, score = c("Identity", "Ranks"), only.mean = FALSE,
  plot = TRUE, FAP = 0.05, seed = 11642257, L = 10000, limits = NA)

mchangepoint(x, subset, score = c("Identity", "Signed Ranks", "Spatial Signs",
  "Spatial Ranks", "Marginal Ranks"), only.mean = FALSE,
  plot = TRUE, FAP = 0.05, seed = 11642257, L = 10000, limits = NA) 

changepoint.normal.limits(n, m, score = c("Identity", "Ranks"),
  only.mean = FALSE, FAP = 0.05, seed = 11642257, L = 100000)

mchangepoint.normal.limits(p, n, m, score = c("Identity", "Signed Ranks", "Spatial Signs",
  "Spatial Ranks", "Marginal Ranks"), only.mean = FALSE,
  FAP = 0.05, seed = 11642257, L = 100000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    \code{changepoint}: a nxm numeric matrix or a numeric vector of
    length m.

    \code{mchangepoint}: a pxnxm data numeric array
    or a pxm numeric vector.

    See below, for the meaning of p, n and m.
  }
  \item{p}{
    integer: number of monitored variables.
  }
  \item{n}{
    integer: size of each subgroup (number of observations gathered at
    each time point).
  }
  \item{m}{
    integer: number of subgroups (time points).
  }
  \item{subset}{
    an optional vector specifying a subset of subgroups/time points to be used
  }
  \item{score}{
    character: the transformation to use; see \code{\link{mshewhart}}.
  }
  \item{only.mean}{
    logical; if \code{TRUE} only a location change-point is searched.
  }
  \item{plot}{
    logical; if \code{TRUE}, the control statistic is displayed. 
  }
  \item{FAP}{
    numeric (between 0 and 1): the desired false alarm probability.
  }
  \item{seed}{
    positive integer; if not \code{NA}, the RNG's state is resetted
    using \code{seed}. The current \code{.Random.seed} will be
    preserved.     Unused by \code{mshewhart} when \code{limits} is not \code{NA}.
  }
  \item{L}{
    positive integer: the number of Monte Carlo replications used to
    compute the control limits.  Unused by \code{changepoint} and
    \code{mchangepoint} when \code{limits} is not \code{NA}.
  }
  \item{limits}{
    numeric: a precomputed vector of length m containing the control limits.
  }
}
\details{
  After an optional rank transformation (argument \code{score}),
  \code{changepoint} and \code{mchangepoint} compute,
  for \eqn{\tau=2,\ldots,m}, the normal likelihood ratio test statistics
  for verifying whether the mean and dispersion (or only the mean when
  \code{only.mean=TRUE}) are the same before and after \eqn{\tau}.
  See Sullivan and Woodall (1999, 2000) and Qiu (2013), Chapter 6 and
  Section 7.5. 

  Note that 
  the control statistic is equivalent to that proposed by
  Lung-Yut-Fong \emph{et al.} (2011)
  when \code{score="Marginal Ranks"} and \code{only.mean=TRUE}.

  As suggested by Sullivan and Woodall (1999, 2000), 
  control limits proportional to the
  in-control mean of the likelihood ratio test statistics
  are used. Further, when \code{plot=TRUE}, the control
  statistics divided by the time-varying control limits
  are plotted with a \dQuote{pseudo-limit} equal to one.

  When \code{only.mean=FALSE}, the decomposition of the
  likelihood ratio test statistic suggested
  by Sullivan and Woodall (1999, 2000)
  for diagnostic purposes is also
  computed, and optionally plotted.
}
\value{
  \code{changepoint} and \code{mchangepoint} return an
  invisible list with elements
  \item{\code{glr} }{
    control statistics.
  }
  \item{\code{mean}, \code{dispersion} }{decomposition
  of the control statistics in the two parts due to changes in the mean and
  dispersion, respectively. These elements are present only when
  \code{only.mean=FALSE}.}
  \item{\code{limits} }{control limits.}
  \item{\code{score}, \code{only.mean}, \code{FAP},
    \code{L}, \code{seed} }{ input arguments.}

  \code{changepoint.normal.limits} and \code{mchangepoint.normal.limits} return a numeric vector
  containing the control limits.
}
\note{
  \enumerate{
    \item
    When \code{limits} is \code{NA}, \code{changepoint}
    and \code{mchangepoint}
    compute the control limits by permutation.
    The resulting control charts are distribution-free.
    \item
    Pre-computed limits, like those computed using
    \code{changepoint.normal.limits} and
    \code{mchangepoint.normal.limits},  
    are recommended \emph{only} for univariate data when \code{score=Ranks}.
    Indeed, in all the other cases, the resulting control
    chart will not be distribution-free.
    \item
    However, note that, when \code{score} is \code{Signed Ranks}, \code{Spatial
    Signs}, \code{Spatial Ranks} 
    the normal-based control limits are distribution-free in the class
    of all multivariate elliptical distributions.
  }
}
\author{
\packageAuthor{dfphase1}.
}
\references{
  A. Lung-Yut-Fong, C. Lévy-Leduc, O. Cappé O (2011) \dQuote{Homogeneity and change-point detection tests
    for multivariate data using rank statistics}. arXiv preprint arXiv:11071971.
  
  P. Qiu (2013) \emph{Introduction to Statistical Process Control}. Chapman & Hall/CRC Press.
  
  J. H. Sullivan, W. H. Woodall (1996) \dQuote{A control chart for preliminary analysis of individual observa-
  tions}. \emph{Journal of Quality Technology}, \bold{28}, pp. 265--278.
  
  J. H. Sullivan, W. H. Woodall (2000) \dQuote{Change-point detection of mean vector or covariance matrix
  shifts using multivariate individual observations}. \emph{IIE Transactions}
  \bold{32}, pp. 537--549.

}
\examples{
data(gravel)
changepoint(gravel[1,,])
mchangepoint(gravel)
mchangepoint(gravel,score="Signed Ranks")
}
\keyword{htest}
\keyword{multivariate}
\keyword{statistical process monitoring}
\keyword{statistical process control}
\keyword{distribution-free methods}
