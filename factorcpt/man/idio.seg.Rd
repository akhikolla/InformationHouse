\name{idio.seg}
\alias{idio.seg}

\title{
Multiple change-point detection for the idiosyncratic components
}

\description{
First generating the panel of statistics via wavelet-based filtering of the estimated
idiosyncratic components, it applies the Double CUSUM Binary Segmentation in combination
with the bootstrap generated thresholds
to estimate the multiple change-points in the idiosyncratic components.
}

\usage{
idio.seg(gfm, q, scales = NULL, diag = TRUE, sig.lev = 0.05, rule = NULL, B = 200,
        p = NULL, dw = NULL, mby = NULL, tby = NULL, do.parallel = TRUE)
}

\arguments{
  \item{gfm}{
  a \code{\link{get.factor.model}} object with estimates of the factor structure
  }
  \item{q}{
  the number of factors
  }
  \item{scales}{
  see \code{scales} in \code{\link{factor.seg.alg}}
  }
  \item{diag}{
  see \code{idio.diag} in \code{\link{factor.seg.alg}}
  }
  \item{sig.lev}{
  see \code{sig.lev} in \code{\link{factor.seg.alg}}
  }
  \item{rule}{
  the height of a binary tree for change-point analysis,
  see the Appendix of Barigozzi, Cho & Fryzlewicz (2016)
  }
  \item{B}{
  the size of bootstrap samples
  }
  \item{p}{
  see \code{p} in \code{\link{factor.seg.alg}}
  }
  \item{dw}{
  see \code{dw} in \code{\link{factor.seg.alg}}
  }
  \item{mby}{
  see \code{dmby} in \code{\link{func_dc_by}}
  }
  \item{tby}{
  see \code{dtby} in \code{\link{func_dc_by}}
  }
  \item{do.parallel}{
  see \code{do.parallel} in \code{\link{factor.seg.alg}}
  }
}

\value{
  \item{tree}{
  a list containing the binary tree grown for change-point analysis
  on the idiosyncratic components; each node contains its index, the start and
  end of the segment as well as the estimated change-point, the test statistic
  and the proportion of bootstrap statistics smaller than the test statistic.
  }
  \item{est.cps}{
  a vector of change-points estimated for the idiosyncratic components;
  adjusted for the bias introduced by the wavelet filtering
  }
  \item{p}{
  the reciprocal of the average block size selected for the estimated idiosyncratic components
  }
}

\references{
M. Barigozzi, H. Cho and P. Fryzlewicz (2016) Simultaneous multiple change-point and factor analysis for high-dimensional time series, Preprint. \cr

H. Cho (2016) Change-point detection in panel data via double CUSUM statistic. Electronic Journal of Statistics. 10: 2000-2038.
}

\author{
Haeran Cho
}

\keyword{internal}
