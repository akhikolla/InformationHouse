\name{factor.seg.alg}
\alias{factor.seg.alg}

\title{
Simultaneous multiple change-point and factor analysis for high-dimensional time series
}

\description{
The function applies the two-stage methodology for multiple change-point detection under factor modelling. It first transforms the input time series into panels of statistics that contain the change-points in the second-order structure of the common and idiosyncratic components, as change-points in their `means', to which the Double CUSUM Binary Segmentation algorithm is applied in the second stage. The function returns change-point estimates from the common and idiosyncratic components separately.
}

\usage{
factor.seg.alg(x, r=NULL, bn.op=2, sig.lev=.05, max.q=NULL, q.seq=NULL,
              qlen=5, qby=0, dw=NULL, p=NULL, B=200, scales=NULL,
              rule=NULL, mby=NULL, tby=NULL, idio.diag=FALSE,
              do.parallel=TRUE, no.proc=2)
}

\arguments{
  \item{x}{
  input time series matrix, with each row representing a time series
  }
  \item{r}{
  the number of factors, if \code{r=NULL}, screening over a range of factor number candidates
  is performed as described in the paper
  }
  \item{bn.op}{
  an index number for the information criterion-based estimator of Bai and Ng (2002)
  for the number of factors, the default value \code{bn.op=2} is used in Barigozzi, Cho & Fryzlewicz (2016)
  }
  \item{sig.lev}{
  sets the level of significance for drawing the bootstrap-based threshold
  }
  \item{max.q}{
  the maximum number of factors, if \code{max.q=NULL}, a default value is chosen as described in the paper
  }
  \item{q.seq}{
  a vector of factor number candidates; if \code{q.seq=NULL}, it is chosen as described in
  Barigozzi, Cho & Fryzlewicz (2016) using \code{qlen} or \code{qby}
  }
  \item{qlen}{
  specifies the length of the sequence of factor number candidates
  }
  \item{qby}{
  specifies the increment of the sequence of factor number candidates when \code{qlen=0};
  if both \code{qlen=0} and \code{qby=0}, the sequence of factor number candidates is generated
  as if \code{qby=1}
  }
  \item{dw}{
  trims off the interval of consideration in the binary segmentation algorithm
  and determines the minimum length of a stationary segment;
  if \code{dw=NULL}, a default value is chosen as described in the Appendix of
  Barigozzi, Cho & Fryzlewicz (2016)
  }
  \item{p}{
  inverse of the average length of blocks in the stationary bootstrap procedure,
  if \code{p=NULL}, a data-driven value is chosen as described in the Appendix of
  Barigozzi, Cho & Fryzlewicz (2016)
  }
  \item{B}{
  the number of bootstrap samples for threshold selection
  }
  \item{scales}{
  the number of wavelet scales used for wavelet filtering of
  the common and idiosyncratic components estimated via PCA
  }
  \item{rule}{
  the depth of a binary tree for change-point analysis, see the
  Appendix of Barigozzi, Cho & Fryzlewicz (2016)
  }
  \item{mby}{
  see \code{dmby} in \code{\link{func_dc_by}}
  }
  \item{tby}{
  see \code{dtby} in \code{\link{func_dc_by}}
  }
  \item{idio.diag}{
  if \code{idio.diag=TRUE}, only the diagonal wavelet-transform is employed
  in order to generate the panel of statistics from the idiosyncratic components
  }
  \item{do.parallel}{
  if \code{do.parallel=TRUE}, a set of copies of R running in parallel are created
  and used for bootstrap procedure
  }
  \item{no.proc}{
  sets the number of processes to be spawned when \code{do.parallel=TRUE},
  }
}

\value{
    \item{cs.list}{a list of objects returned by an internal function \code{\link{common.seg}},
    which contain the results from performing change-point analysis on the common components
    estimated with a range of factor number candidates}
    \item{r}{the factor number selected from performing change-point analysis on the common component}
    \item{common.est.cps}{change-points detected from the common component estimated with \code{r} as the factor number}
    \item{idio.seg.res}{an object returned by an internal function \code{\link{idio.seg}},
    which contains the results from performing change-point analysis on the idiosyncratic
    component estimated with \code{r} as the factor number}
    \item{idio.est.cps}{change-points detected from the idiosyncratic component}
    \item{gfm}{factor structure of \code{x} estimated from \code{\link{get.factor.model}}}
    \item{q.seq}{a vector containing the range of factor number candidates}
}

\references{
J. Bai and S. Ng (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221. \cr

M. Barigozzi, H. Cho and P. Fryzlewicz (2016) Simultaneous multiple change-point and factor analysis for high-dimensional time series, Preprint. \cr

H. Cho (2016) Change-point detection in panel data via double CUSUM statistic. Electronic Journal of Statistics. 10: 2000-2038.
}

\author{ Haeran Cho }

\examples{
n <- 50; T <- 200
e <- matrix(rnorm(n*T), nrow=n) # idiosyncratic components
r <- 3 # factor number
Lam <- matrix(rnorm(n*r, 1, 1), nrow=n) # loadings
f <- matrix(rnorm(r*T), nrow=r) # factors
chi <- e*0 # common component
chp <- T/2 # change-point
chi[, 1:chp] <- Lam\%*\%f[, 1:chp]
Lam <- Lam + matrix(rnorm(n*r, 0, sqrt(2)), nrow=n) # new loadings
chi[, (chp+1):T] <- Lam\%*\%f[, (chp+1):T]
x <- chi + sqrt(r)*e

fsa <- factor.seg.alg(x, idio.diag=TRUE)
fsa$common.est.cps
fsa$idio.est.cps
fsa$q.seq
}
