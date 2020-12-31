#' Efficient calculation of covariance matrices for Raster* objects
#'
#' \code{parCov} efficiently calculates the covariance of Raster* objects,
#' taking advantage of parallel processing and pulling data into memory only as
#' necessary. For large datasets with lots of variables, calculating the covariance
#' matrix rapidly becomes unwieldy, as the number of calculations required grows
#' quadratically with the number of variables.
#'
#' @param x Raster* object, typically a brick or stack
#' @param y NULL (default) or a Raster* object with the same extent and resolution
#'   as \code{x}
#' @param ... additional arguments, including any of the following:
#' @param w optional Raster* object of weights for a weighted covariance matrix
#' @param sample logical. If \code{TRUE}, the sample covariance is calculated
#'   with a denominator of $n-1$
#' @param progress logical. If \code{TRUE}, messages and progress bar will be
#'   printed
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Number of CPU cores to utilize for parallel processing
#' @param cl optional cluster object
#' @param keep.open logical. If \code{TRUE} and \code{parallel = TRUE}, the
#'   cluster object will not be closed after the function has finished
#'
#' @examples
#' mat1 <- parCov(climdat.hist)
#'
#' # correlation matrix
#' Z <- parScale(climdat.hist)
#' mat2 <- parCov(Z)
#'
#' # covariance between two Raster* objects
#' mat3 <- parCov(x = climdat.hist, y = climdat.fut)
#'
#' @return Returns a matrix with the same row and column names as the layers of
#'   \code{x}. If \code{y} is supplied, then the covariances between the layers
#'   of \code{x} and the layers of code{y} are computed.
#'
#' @details This function is designed to work similarly to the
#'   \code{\link[stats]{cov}} and the \code{\link[raster]{layerStats}}
#'   functions, with two major differences. First, \code{parCov} allows you to
#'   calculate the covariance between two different Raster* objects, whereas
#'   \code{layerStats} does not. Second, \code{parCov} can (optionally) compute
#'   each element of the covariance matrix in parallel, offering a dramatic
#'   improvement in computation time for large Raster* objects.
#'
#'   The raster layer of weights \code{w} should contain raw weights as values,
#'   and should \emph{not} be normalized so that \code{sum(w) = 1}. This is
#'   necessary for computing the sample covariance, whose formula contains
#'   \code{sum(w) - 1} in its denominator.
#'
#' @seealso \code{\link[stats]{cov}}, \code{\link[raster]{layerStats}}
#'
#' @export
#' @importFrom pbapply pbsapply pboptions
#' @importFrom foreach '%dopar%'
#' @importFrom doSNOW registerDoSNOW
#' @importFrom snow makeCluster clusterExport
#' @importFrom parallel detectCores

setGeneric("parCov", function(x, y, ...){
  standardGeneric("parCov")})

#' @rdname parCov
setMethod("parCov",
          signature(x = "Raster", y = "missing"),
          function(x, w = NULL, sample = TRUE, progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE){

            if (canProcessInMemory(x) && !parallel) {
              dat <- values(x)
              mat <- stats::cov(dat, method = "pearson", use = "pairwise.complete.obs")
              return(mat)
            }
            nl <- nlayers(x)
            mat <- matrix(NA, nrow = nl, ncol = nl)
            colnames(mat) <- rownames(mat) <- names(x)

            ii <- rep(1, nl)
            for(i in 2:nl) ii <- c(ii, rep(i, each = (nl - i + 1)))
            jj <- 1:nl
            for(i in 2:nl) jj <- c(jj, i:nl)
            s <- 1:length(ii)

            if (!parallel || n == 1) {
              if (progress) {
                pboptions(char = "-", txt.width = NA, type = "txt")
                result <- pbsapply(s, function(p) do.call(.covij, list(x = subset(x, ii[p]), y = subset(x, jj[p]), w = w, sample = sample)))
              } else if (!progress) {
                result <- sapply(s, function(p) do.call(.covij, list(x = subset(x, ii[p]), y = subset(x, jj[p]), w = w, sample = sample)))
              }
            }

            if (parallel && n > 1) {
              if (!keep.open) on.exit(closeAllConnections())
              if(!is.numeric(n) && is.null(cl)) {
                n <- min(detectCores() - 1, floor(length(s)/2))
                if (progress) message('incorrect number of cores specified, using ', n)
              } else if(is.null(cl) && n > detectCores()) {
                n <- min(detectCores() - 1, floor(length(s)/2))
                if (progress) message('too many cores specified, using ', n)
              }
              w <- w
              if (is.null(cl)) cl <- makeCluster(getOption("cl.cores", n))
              clusterExport(cl, c(".covij", "raster", "cellStats", "x", "ii", "jj", "s", "w", "canProcessInMemory", "values", "sample", "subset"),
                            envir = environment())
              registerDoSNOW(cl)
              if (progress) {
                pb <- txtProgressBar(min = 0, max = length(s), style = 3, char = "-")
                progress <- function(n) setTxtProgressBar(pb, n)
                opts <- list(progress = progress)
                result <- foreach::foreach(p = s, .combine=c, .options.snow = opts) %dopar% {
                  do.call(.covij, list(x = subset(x, ii[p]), y = subset(x, jj[p]), w = w, sample = sample))
                }
                close(pb)
              } else if(!progress) {
                result <- foreach::foreach(p = s, .combine=c) %dopar% {
                  do.call(.covij, list(x = subset(x, ii[p]), y = subset(x, jj[p]), w = w, sample = sample))
                }
              }
              if (!keep.open || is.null(cl)) snow::stopCluster(cl)
            }

            for (p in s) {
              mat[ii[p], jj[p]] <- mat[jj[p], ii[p]] <- result[p]
            }

            return(mat)
          }
)

#' @rdname parCov
setMethod("parCov",
          signature(x = "Raster", y = "Raster"),
          function(x, y, w = NULL, sample = TRUE, progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE){

            if (canProcessInMemory(x) & !parallel) {
              x.dat <- values(x)
              y.dat <- values(y)
              mat <- stats::cov(x.dat, y.dat, method = "pearson", use = "pairwise.complete.obs")
              return(mat)
            }

            nlx <- nlayers(x)
            nly <- nlayers(y)
            mat <- matrix(NA, nrow = nlx, ncol = nly)
            rownames(mat) <- names(x)
            colnames(mat) <- names(y)
            z <- expand.grid(1:nlx, 1:nly)
            s <- 1:nrow(z)

            if (!parallel | n == 1) {
              if (progress) {
                pboptions(char = "-", txt.width = NA, type = "txt")
                result <- pbsapply(s, function(p) do.call(.covij, list(x = subset(x, z[p, 1]), y = subset(y, z[p, 2]), w = w, sample = sample)))
              } else if (!progress){
                result <- sapply(s, function(p) do.call(.covij, list(x = subset(x, z[p, 1]), y = subset(y, z[p, 2]), w = w, sample = sample)))
              }
            }

            if (parallel && n > 1) {
              if (!keep.open) on.exit(closeAllConnections())
              if(!is.numeric(n) && is.null(cl)) {
                n <- min(detectCores() - 1, floor(length(s)/2))
                if (progress) message('incorrect number of cores specified, using ', n)
              } else if(is.null(cl) && n > parallel::detectCores()) {
                n <- min(detectCores() - 1, floor(length(s)/2))
                if (progress) message('too many cores specified, using ', n)
              }
              w <- w
              if (is.null(cl)) cl <- snow::makeCluster(getOption("cl.cores", n))
              snow::clusterExport(cl, c(".covij", "raster", "cellStats", "x", "y", "z", "s", "w", "canProcessInMemory", "values", "sample", "subset"),
                                  envir = environment())
              doSNOW::registerDoSNOW(cl)
              if (progress) {
                pb <- txtProgressBar(min = 0, max = length(s), style = 3, char = "-")
                progress <- function(n) setTxtProgressBar(pb, n)
                opts <- list(progress = progress)
                result <- foreach::foreach(p = s, .combine=c, .options.snow = opts) %dopar% {
                  do.call(.covij, list(x = subset(x, z[p, 1]), y = subset(y, z[p, 2]), w = w, sample = sample))
                }
                close(pb)
              } else if (!progress) {
                result <- foreach::foreach(p = s, .combine=c) %dopar% {
                  do.call(.covij, list(x = subset(x, z[p, 1]), y = subset(y, z[p, 2]), w = w, sample = sample))
                }
              }
              if (!keep.open || is.null(cl)) snow::stopCluster(cl)
            }

            for (p in s) {
              mat[ z[p,2], z[p,1] ] <- result[p]
              #if(nlx > 1 & nly > 1) mat[ z[p,1], z[p,2] ] <- mat[ z[p,2], z[p,1] ]
            }
            return(mat)
          }
)

#' @keywords internal
.expand.grid.unique <- function(x, y){

  nx <- length(x)
  ny <- length(y)

  if(ny == 1){
    dat <- cbind(x, y)
  } else if(nx <= ny) {
    dat <- NULL
    for(i in 1:nx){
      for(j in i:ny){
        dat <- rbind(dat, c(x[i], y[j]))
      }
    }
  } else if(nx > ny) {
    dat <- NULL
    for(j in 1:ny){
      for(i in j:nx){
        dat <- rbind(dat, c(x[i], y[j]))
      }
    }
  }
  return(dat)
}
