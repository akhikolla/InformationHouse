#' Multiscale MOSUM algorithm with bottom-up merging
#' 
#' Multiscale MOSUM procedure with symmetric bandwidths combined with
#' bottom-up bandwidth-based merging.
#' @param x input data (a \code{numeric} vector or an object of classes \code{ts} and \code{timeSeries})
#' @param G a vector of (symmetric) bandwidths, given as either integers less than \code{length(x)/2},
#' or numbers between \code{0} and \code{0.5} describing the moving sum bandwidths relative to \code{length(x)}.
#' If the smallest bandwidth is smaller than \code{min(20, 0.05*length(x))} 
#' (\code{0.05} if relative bandwidths are given) and \code{threshold = "critical.value"},
#' it generates a warning message
#' @param threshold string indicating which threshold should be used to determine significance.
#' By default, it is chosen from the asymptotic distribution at the given significance level \code{alpha}.
#' Alternatively, it is possible to parse a user-defined function with \code{threshold.function}
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff \code{threshold = "critical.value"}
#' @param threshold.function function object of form \code{function(G, length(x), alpha)}, to compute a
#' threshold of significance for different bandwidths G; use iff \code{threshold = "custom"}
#' @param eta see \link[mosum]{mosum}
#' @param do.confint flag indicating whether to compute the confidence intervals for change-points
#' @param level use iff \code{do.confint = TRUE}; a numeric value (\code{0 <= level <= 1}) with which
#' \code{100(1-level)\%} confidence interval is generated
#' @param N_reps use iff \code{do.confint = TRUE}; number of bootstrap replicates to be generated
#' @param ... further arguments to be passed to the \link[mosum]{mosum} calls
#' @return S3 object of class \code{multiscale.cpts}, which contains the following fields:
#'    \item{x}{input data}
#'    \item{cpts}{estimated change-points}
#'    \item{cpts.info}{data frame containing information about estimated change-points}
#'    \item{pooled.cpts}{set of change-point candidates that have been considered by the algorithm}
#'    \item{G}{bandwidths}
#'    \item{threshold, alpha, threshold.function}{input parameters}
#'    \item{eta}{input parameters}
#'    \item{do.confint}{input parameter}
#'    \item{ci}{object of class \code{cpts.ci} containing confidence intervals for change-points iff \code{do.confint = TRUE}}
#' @details See Algorithm 1 in the first referenced paper for a comprehensive
#' description of the procedure and further details.
#' @references A. Meier, C. Kirch and H. Cho (2019)
#' mosum: A Package for Moving Sums in Change-point Analysis. \emph{To appear in the Journal of Statistical Software}.
#' @references M. Messer et al. (2014)
#' A multiple filter test for the detection of rate changes in renewal processes with varying variance.
#' \emph{The Annals of Applied Statistics}, Volume 8, Number 4, pp. 2027-2067.
#' @examples 
#' x1 <- testData(lengths = c(100, 200, 300, 300), 
#' means = c(0, 1, 2, 2.7), sds = rep(1, 4), seed = 123)$x
#' mbu1 <- multiscale.bottomUp(x1)
#' plot(mbu1)
#' summary(mbu1)
#' 
#' x2 <- testData(model = "mix", seed = 1234)$x
#' threshold.custom <- function(G, n, alpha) {
#' mosum.criticalValue(n, G, G, alpha) * log(n/G)^0.1
#' }
#' mbu2 <- multiscale.bottomUp(x2, G = 10:40, threshold = "custom",
#' threshold.function = threshold.custom)
#' plot(mbu2)
#' summary(mbu2)
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib mosum, .registration = TRUE
#' @export
multiscale.bottomUp <- function(x, G=bandwidths.default(length(x), G.min = max(20, ceiling(0.05*length(x)))), 
                                threshold = c('critical.value', 'custom')[1], 
                                alpha=0.1, threshold.function = NULL, eta=0.4, 
                                do.confint=FALSE, level=0.05, N_reps=1000, ...) {
  n <- length(x)

  if (class(G) == 'integer' || class(G) == 'numeric') {
    grid <- multiscale.grid(G, method='concatenate')
   } else if (class(G) == 'multiscale.grid'){
    if (any(apply(G$grid, 1, diff) != 0)) {
      stop("Expecting a grid of symmetric bandwidths")
    }
    grid <- G
  } else stop('Expecting a vector of numbers')
  abs.bandwidth <- all(grid$grid>=1)
  
  if (abs.bandwidth) {
    GRID_THRESH <- max(20, 0.05*n)
  } else {
    GRID_THRESH <- 0.05
  }
  
  if (threshold == 'critical.value' & min(grid$grid) < GRID_THRESH) {
     warning('Smallest bandwidth in grid is relatively small (in comparison to n), \n increase the smallest bandwidth or use multiscale.localPrune instead')
  }
  if (threshold != 'critical.value' && threshold != 'custom') {
    stop('threshold must be either \'critical.value\' or \'custom\'')
  }

  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(eta <= 1 & eta > 0) 
  stopifnot(!do.confint || N_reps>0)
  
  # Retreive change-point candidates from all bandwidths.
  cpts.complete <- numeric(0)
  bandwidths.complete <- integer(0)
  pValues.complete <- numeric(0)
  jumps.complete <- numeric(0)

  for (i in seq_len(nrow(grid$grid))) {
    G <- grid$grid[[i, 1]]
    if (threshold == 'critical.value') {
      m <- mosum(x, G, ..., threshold='critical.value', alpha=alpha, criterion='eta', eta=eta)
    } else {
      threshold_val <- threshold.function(G, n, alpha)
      m <- mosum(x, G, ... , threshold='custom', threshold.custom=threshold_val, alpha=alpha, criterion='eta', eta=eta)
    }
    if(!abs.bandwidth) G <- floor(G*n)
    cpts <- m$cpts
    cpts.complete <- c(cpts.complete, cpts)
    bandwidths.complete <- c(bandwidths.complete, rep(G, length(cpts)))
    pValues.complete <- c(pValues.complete, mosum.pValue(m$stat[cpts], n, G))
    jumps.complete <- c(jumps.complete, m$stat[cpts]*sqrt(2/G))
  }

  # Merge candidates.
  points <- numeric(0)
  bandwidths <- numeric(0)
  pValues <- numeric(0)
  jumps <- numeric(0)
  cptsInOrder <- seq_len(length(cpts.complete))
  for (i in cptsInOrder) {
    p <- cpts.complete[[i]]
    G <- bandwidths.complete[[i]]
    pVal <- pValues.complete[[i]]
    jmp <- jumps.complete[[i]]
    if (suppressWarnings(min(abs(p-points))) >= eta*G) { # Note: min(empty_list) = Inf
      points <- c(points, p)
      bandwidths <- c(bandwidths, G)
      pValues <- c(pValues, pVal)
      jumps <- c(jumps, jmp)
    }
  }
  cpts.merged <- data.frame(cpts = points, G.left = bandwidths, G.right = bandwidths, 
                            p.value = pValues, jump = jumps)
  cpts <- cpts.merged[order(cpts.merged$cpts), ]
  G <- as.vector(grid$grid[, 1])
  if(!abs.bandwidth) G <- floor(n*G)
  
  ret <- structure(list(x=x,
                        cpts=as.numeric(cpts[, 1]), 
                        cpts.info=cpts,
                        pooled.cpts=sort(unique(cpts.complete)), 
                        G=G,
                        alpha=alpha,
                        threshold=threshold,
                        threshold.function=threshold.function,
                        criterion='eta',
                        eta=eta,
                        do.confint=FALSE,
                        ci=NA), # note 
                   class='multiscale.cpts')
  if (do.confint) {
    ret$ci <- confint.multiscale.cpts(ret, level=level, N_reps=N_reps)
    ret$do.confint <- TRUE
  }
  ret
}
