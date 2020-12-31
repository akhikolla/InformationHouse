#' Brier Score decomposition
#'
#' @description Return decomposition of the Brier Score into Reliability, Resolution and Uncertainty, and estimated standard deviations
#' @aliases BrierScoreDecomposition BrierDecomp
#' @param p vector of forecast probabilities
#' @param y binary observations, y[t]=1 if an event happens at time t, and y[t]=0 otherwise
#' @param bins binning to estimate the calibration function (see Details), default: 10
#' @param bias.corrected logical, default=FALSE, whether the standard (biased) decomposition of Murphy (1973) should be used, or the bias-corrected decomposition of Ferro (2012)
#' @return Estimators of the three components and their estimated standard deviations are returned as a 2*3 matrix.
#' @details
#' To estimate the calibration curve, the unit line is categorised into discrete bins, provided by the `bins` argument. If `bins` is a single number, it specifies the number of equidistant bins. If `bins` is a vector of values between zero and one, these values are used as the bin-breaks. 
#' @examples
#' data(eurotempforecast)
#' BrierDecomp(rowMeans(ens.bin), obs.bin, bins=3, bias.corrected=TRUE)
#' @seealso ReliabilityDiagram
#' @references
#' Murphy (1973): A New Vector Partition of the Probability Score. J. Appl. Met. \doi{10.1175/1520-0450(1973)012<0595:ANVPOT>2.0.CO;2}
#'
#' Ferro and Fricker (2012): A bias-corrected decomposition of the Brier score. QJRMS. \doi{10.1002/qj.1924}
#'
#' Siegert (2013): Variance estimation for Brier Score decomposition. QJRMS. \doi{10.1002/qj.2228}
#' @export

BrierDecomp <- function(p, y, bins=10, bias.corrected=FALSE) {

  n <- length(p)

  # binning, either number of bins or bin breaks are specified
  if (length(bins) == 1) {
    n.bins <- floor(bins)
    p.breaks <- seq(0, 1, length.out=n.bins+1) 
  } else {
    n.bins <- length(bins) - 1
    bins <- sort(bins)
    stopifnot(min(bins)<= 0 & max(bins) >= 1)
    p.breaks <- bins
  }


  #p.binning (vector of length n with entries equal to bin-no. of p_n)
  p.binning <- cut(p, breaks=p.breaks, include.lowest=TRUE, ordered_result=TRUE)
  p.binning <- as.numeric(p.binning)

  ##construct matrices and column sums
  m.ind <- matrix(c(seq(n), p.binning), nrow=n)
  m.A <- matrix(0, nrow=n, ncol=n.bins)
  m.A[m.ind] <- 1
  cs.A <- colSums(m.A)
  m.B <- matrix(0, nrow=n, ncol=n.bins)
  m.B[m.ind] <- y
  cs.B <- colSums(m.B)
  m.C <- matrix(0, nrow=n, ncol=n.bins)
  m.C[m.ind] <- p
  cs.C <- colSums(m.C)
  m.Y <- matrix(y, ncol=1)
  cs.Y <- colSums(m.Y)

  ## sets of indices for which colSums(A) > 0, resp. > 1
  ## in order to avoid division by zero
  d0 <- which(cs.A < 1)
  d1 <- which(cs.A < 2)
  d0.set <- function(x) {
    if (length(d0) > 0) {
      x <- x[-d0]
    }
    x
  }
  d1.set <- function(x) {
    if (length(d1) > 0) {
      x <- x[-d1]
    }
    x
  }

  ## calculate estimators
  rel <- 1/n * sum( d0.set((cs.B - cs.C)^2 / cs.A) )
  res <- 1/n * sum( d0.set(cs.A * (cs.B / cs.A - cs.Y / n)^2) )
  unc <- cs.Y * (n - cs.Y) / n / n

  ## correction factors for bias correction
  corr.s <- 1 / n * sum(d1.set(cs.B * (cs.A - cs.B) / (cs.A * (cs.A - 1))))
  corr.t <- cs.Y * (n - cs.Y) / n / n / (n-1)

  # avoid rel<0, res<0, res>1, and unc>1/4
  alpha <- min(c(rel/corr.s, 
                 max(c(res / (corr.s - corr.t), 
                       (res - 1) / (corr.s - corr.t))),
                 (1 - 4 * unc) / (4 * corr.t),
                 1)
              )
  if (!is.finite(alpha)) {
    alpha <- 0
  }

  rel2 <- rel - alpha * corr.s
  res2 <- res - alpha * corr.s + alpha * corr.t
  unc2 <- unc + alpha * corr.t
                  


  ###################################################
  # variances of reliability, resolution, uncertainty
  # estimated by propagation of error
  ###################################################

  #####
  # REL
  #####
  d.rel.d.a <- -1 * (cs.B - cs.C) * (cs.B - cs.C) / (n * cs.A * cs.A)
  d.rel.d.b <- 2 * (cs.B - cs.C) / (n * cs.A)
  d.rel.d.c <- -2 * (cs.B - cs.C) / (n * cs.A)
  if (length(d0) > 0) {
    d.rel.d.a[d0] <- 0
    d.rel.d.b[d0] <- 0
    d.rel.d.c[d0] <- 0
  }
  jacobian.rel <- matrix(c(d.rel.d.a, d.rel.d.b, d.rel.d.c), nrow=1)

  m.X <- cbind(m.A, m.B, m.C)
  cov.X.rel <- crossprod(scale(m.X, center=TRUE, scale=FALSE))

  var.rel <- jacobian.rel %*% cov.X.rel %*% t(jacobian.rel)

  #####
  # RES
  #####
  d.res.d.a <- -1 / n * (cs.B / cs.A - cs.Y / n) * 
                      (cs.B / cs.A + cs.Y / n)
  d.res.d.b <- 2 / n * (cs.B / cs.A - cs.Y / n)
  d.res.d.y <- 0
  if (length(d0) > 0) {
    d.res.d.a[d0] <- 0
    d.res.d.b[d0] <- 0
  }
  jacobian.res <- matrix(c(d.res.d.a, d.res.d.b, d.res.d.y), nrow=1)

  m.X <- cbind(m.A, m.B, m.Y)
  cov.X.res <- crossprod(scale(m.X, center=TRUE, scale=FALSE))
  
  var.res <- jacobian.res %*% cov.X.res %*% t(jacobian.res)

  #####
  # UNC
  #####
  var.unc <- (1 - 2 * cs.Y / n) * (1 - 2 * cs.Y / n) / n / n * 
             crossprod(scale(m.Y, center=TRUE, scale=FALSE))

  ######
  # REL'
  ######
  d.rel2.d.a <- -1 * ((cs.B - cs.C) * (cs.B - cs.C) + 
                      (cs.B * cs.B) / (cs.A - 1) - 
                      cs.A * cs.B * (cs.A - cs.B) / 
		      (cs.A - 1) / (cs.A - 1)
                     ) / (n * cs.A * cs.A)
  d.rel2.d.b <- (2 * cs.B - 1) / (n * cs.A - n) - 2 * cs.C / (n * cs.A)
  d.rel2.d.c <- -2 * (cs.B - cs.C) / (n * cs.A)
  if (length(d1) > 0) {
    d.rel2.d.a[d1] <- 0
    d.rel2.d.b[d1] <- 0
    d.rel2.d.c[d1] <- 0
  }
  jacobian.rel <- matrix(c(d.rel2.d.a, d.rel2.d.b, d.rel2.d.c), nrow=1)

  var.rel2 <- jacobian.rel %*% cov.X.rel %*% t(jacobian.rel)


  ######
  # RES'
  ######

  d.res2.d.a <- -1 / n * (cs.B / cs.A - cs.Y / n) * (cs.B / cs.A + cs.Y / n) +
                 cs.B / (n * cs.A * cs.A * (cs.A - 1) * (cs.A - 1)) *
                 ((cs.A - cs.B) * (cs.A - cs.B) - cs.B * (cs.B - 1))
  d.res2.d.b <- 2 / n * (cs.B / cs.A - cs.Y / n) - (cs.A - 2 * cs.B) /
                (n * cs.A * (cs.A - 1))
  d.res2.d.y <- (n - 2 * cs.Y) / n / n / (n - 1)
  if (length(d1) > 0) {
    d.res2.d.a[d1] <- 0
    d.res2.d.b[d1] <- 0
  }

  jacobian.res2 <- matrix(c(d.res2.d.a, d.res2.d.b, d.res2.d.y), nrow=1)

  var.res2 <- jacobian.res2 %*% cov.X.res %*% t(jacobian.res2)

  ######
  # UNC'
  ######
 
  var.unc2 <- (n - 2 * cs.Y) * (n - 2 * cs.Y) / n / n / (n - 1) / (n - 1) *
              crossprod(scale(m.Y, center=TRUE, scale=FALSE))
  

  #store
  sqrt0 <- function(x) {
    if (x <= 0) {
      return(0)
    } else {
      return(sqrt(x))
    }
  }
  if (bias.corrected == FALSE) {
    ans <- rbind(
      component = c(REL=rel, RES=res, UNC=unc),
      component.sd = c(REL=sqrt0(drop(var.rel)), 
                RES=sqrt0(drop(var.res)), 
                UNC=sqrt0(drop(var.unc)))
    )
  } else {
    ans <- rbind(
      component = c(REL=rel2, RES=res2, UNC=unc2),
      component.sd = c(REL=sqrt0(drop(var.rel2)), 
                RES=sqrt0(drop(var.res2)), 
                UNC=sqrt0(drop(var.unc2)))
    )

  }

  #return
  return(ans)
}

#' @export
BrierScoreDecomposition <- BrierDecomp
