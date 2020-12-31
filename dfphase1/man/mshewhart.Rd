\name{mshewhart}
\alias{mshewhart}
\alias{mshewhart.normal.limits}
\title{
  Multivariate Shewhart-type control charts
}
\description{
  \code{mshewhart} computes, and, optionally,  plots, 
  several Shewhart-type Phase I control charts for detecting
  location and scale changes in multivariate subgrouped data.
  

  \code{mshewhart.normal.limits} pre-computes
  the corresponding control limits when the in-control distribution is
  multivariate normal.
}
\usage{
mshewhart(x, subset, stat = c("T2Var", "T2", "Var", "Depth Ranks"), score = c("Identity",
  "Signed Ranks",  "Spatial Signs", "Spatial Ranks", "Marginal Ranks"),
  loc.scatter = c("Classic", "MCD"), plot = TRUE, FAP = 0.05,
  seed = 11642257, L = 1000, limits = NA)

mshewhart.normal.limits(p, n, m, stat = c("T2Var", "T2", "Var", "Depth Ranks"),
  score = c("Identity", "Signed Ranks",  "Spatial Signs", "Spatial Ranks",
  "Marginal Ranks"), loc.scatter = c("Classic", "MCD"),
  FAP = 0.05, seed = 11642257, L = 100000) 
}
\arguments{
  \item{x}{
    a pxnxm data numeric array
    (n observations gathered at m time points on p variables).
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
  \item{stat}{
    character: control statistic[s] to use; see Details.
  }
  \item{score}{
    character: transformation to use; unused when
  \code{stat=Depth Ranks}; see Details.
  }
  \item{loc.scatter}{
    character: estimates of the multivariate location and scatter
    to use when no preliminary rank transformation is applied.
    Unused when \code{stat} is equal to \code{Depth Ranks} or \code{score} is different
    from \code{Identity}. See Details.
  }
  \item{plot}{
    logical; if \code{TRUE}, control statistic[s] is[are] displayed. 
  }
  \item{FAP}{
    numeric (between 0 and 1): desired false alarm probability.
  }
  \item{seed}{
    positive integer; if not \code{NA}, the RNG's state is resetted
    using \code{seed}. The current \code{.Random.seed} will be
    preserved.     Unused by \code{mshewhart} when \code{limits} is not \code{NA}.
  }
  \item{L}{
    positive integer: number of Monte Carlo replications used to
    compute the control limits.  Unused by \code{mshewhart} when \code{limits} is not \code{NA}.
  }
  \item{limits}{
    numeric: pre-computed vector of control limits.
    This vector should contain \eqn{(A,B)}
    when \code{stat=T2Var}, \eqn{(A)} when \code{stat=T2},
    \eqn{(B)} when \code{stat=Var} and
    \eqn{(C)} when \code{stat=Depth Ranks}.
    See Details for the definition of the critical values
    \eqn{A}, \eqn{B} and \eqn{C}. 
  }
}
\details{
  The implemented control statistics are
  \itemize{
    \item{\code{T2Var}:} combination of the \code{T2} and
      \code{Var} statistics described below.
    \item{\code{T2}:} Hotelling's \eqn{T^2}{T2} control statistics
    (see Montgomery, 2009, equation 11.19, or Qiu, 2013, equation
    7.7) with control limit equal to \eqn{A}.
    \item{\code{Var}:} normal likelihood ratio control statistics
    for detecting changes in the multivariate dispersion
    (see Montgomery, 2009, equation 11.34), with control limit
    equal to \eqn{B}.
    \item{\code{Depth Ranks}:}
    control statistics based on the rank of the Mahalanobis
    depths, proposed by Bell \emph{et. al.}. As suggested Bell \emph{et
      al.}, the Mahalanobis depths are computed using the BACON estimates
    of the multivariate mean vector and the mean of the subgroups sample
    covariance matrices.
    An alarm is signalled if any of the statistics is greater
    than a positive control limit \eqn{C}.
  }
  
  The \code{T2} and \code{Var} control statistics are computed
  \itemize{
    \item{\code{score=Identical:}}
      from the original data standardized
      using either the classical pooled estimates of the mean vector
      and dispersion matrix (Montgomery, 2009, equations 11.14--11.18;
      Qiu, 2013, equations at page 269)   or the highly robust minimum covariance determinant (MCD)
      estimate when argument \code{loc.scatter} is equal to
      \code{Classic} or \code{MCD}, respectively.
    \item{\code{score=Signed Ranks, Spatial Signs, Spatial Ranks,
  Marginal Ranks}:} from a \dQuote{rank} transformation of the original 
  data. In particular, see   Hallin and Paindaveine (2005) for the
  definition of the multivariate signed ranks and Oja (2010) for those of the
  spatial signs, spatial ranks, and marginal ranks. Multivariate signed
  ranks, spatial signs and
  ranks are \dQuote{inner} standardized while marginal ranks are
  \dQuote{outer} standardized (see Oja (2010) for the definition of
  \dQuote{inner} and \dQuote{outer} standardization).

  Note that the \eqn{T^2}{T2} control statistics based on the spatial
  signs corresponds to the control charts suggested by
  Cheng and Shiau (2015).
  }
}
\value{
  \code{mshewhart} returns an invisible list with elements:
  \item{\code{T2} }{\eqn{T^2}{T2} control statistic; this element is present only if
    \code{stat} is \code{T2Var} or \code{T2}.}
  \item{\code{Var} }{\eqn{Var} control statistic; this element is present only if
    \code{stat} is \code{T2Var} or \code{Var}.}
  \item{\code{DepthRanks} }{control statistic based on the rank of
    the Mahalanobis depths; this element is present only if
    \code{stat} is \code{Depth Ranks}.}
  \item{\code{center}, \code{scatter} }{estimates of the
    multivariate location and scatter used to standardized the observations.}
    \item{\code{limits} }{control limits.}
  \item{\code{stat}, \code{score}, \code{loc.scatter}, \code{FAP},
    \code{L}, \code{seed} }{ input arguments.}

  \code{mshewhart.normal.limits} returns a numeric vector
  containing the control limits.
}
\note{
  \enumerate{
    \item
    When \code{limits} is \code{NA}, \code{mshewhart}
    computes the control limits by permutation.
    Then, the resulting control chart is distribution-free.
    \item
    Pre-computed limits, such as  those computed by
    using \code{mshewhart.normal.limits},  are not recommended.
    Indeed, the resulting control
    chart will not be distribution-free.
    \item
    However, when \code{score} is \code{Signed Ranks}, \code{Spatial
      Signs}, \code{Spatial Ranks} or \code{stat} is \code{Depth Ranks},
    the computed control limits are distribution-free in the class
    of all multivariate elliptical distributions.
  }
}
\author{
\packageAuthor{dfphase1}.
}

\references{
  R. C. Bell, L. A. Jones-Farmer, N. Billor (2014) \dQuote{A distribution-free multivariate Phase I location control
chart for subgrouped data from elliptical
  distributions}. \emph{Technometrics}, \bold{56}, pp. 528--538.

  
  C. R. Cheng, J. J. H. Shiau JJH (2015)
  \dQuote{A distribution-free multivariate control chart for Phase I applications}.
    \emph{Quality and Reliability Engineering International}, \bold{31},
    pp. 97--111.

  M. Hallin and D. Paindaveine (2005) \dQuote{Affine-Invariant Aligned
  Rank Tests for the Multivariate General Linear Model with VARMA
  Errors}. \emph{Journal of Multivariate Analysis}, \bold{93}, pp. 122--163.
  
  D. C. Montgomery (2009) \emph{Introduction to Statistical Quality
    Control}, 6th edn. Wiley.

  H. Oja (2010) \emph{Multivariate Nonparametric Methods with R. An Approach Based on Spatial Signs
    and Ranks}. Springer.

  P. Qiu (2013) \emph{Introduction to Statistical Process Control}. Chapman & Hall/CRC Press.
}




\examples{
data(ryan)
mshewhart(ryan)
mshewhart(ryan,subset=-10)
mshewhart(ryan,subset=-c(10,20))
mshewhart(ryan,score="Signed Ranks")
mshewhart(ryan,subset=-10,score="Signed Ranks")
mshewhart(ryan,subset=-c(10,20),score="Signed Ranks")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{statistical process monitoring}
\keyword{statistical process control}
\keyword{distribution-free methods}