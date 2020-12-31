#' Simulate from LSBCLUST Model
#' 
#' Simulate three-way arrays adhering to the LSBCLUST framework (see \code{\link{lsbclust}}).
#' 
#' @param ndata Integer giving the number of data sets to generate with the same underlying
#' parameters.
#' @param nobs Integer giving the number of observations to sample.
#' @param size Vector with two elements giving the number of rows and columns respectively
#' of each simulated observation.
#' @param clustsize A list of length four, with each element containing a vector 
#' of the same length as the corresponding entry in \code{nclust}, indicating the 
#' number of elements to contribute to each sample. Naturally, each of these 
#' vectors must sum to \code{nobs}, or an error will result. Positional matching 
#' are used, in the order "overall", "rows", "columns" and "interactions". If 
#' \code{NULL}, all clusters will be of equal size.
#' @param err_sd The standard deviation of the error distribution, as passed to 
#' \code{\link{rnorm}}
#' @param svmins Vector of minimum values for the singular values 
#' (as passed to \code{\link{simsv}}). Optionally, if all minima are equal,
#' a single numeric value which will be expanded to the correct length.
#' @param svmax The maximum possible singular value (as passed to \code{\link{simsv}})
#' @inheritParams lsbclust
#' @importFrom mvtnorm rmvnorm
#' @export
#' @examples
#' ## Nothing fixed, balanced classes
#' set.seed(1)
#' dat <- rlsbclust(ndata = 1, nobs = 100, size = c(10, 8), nclust = c(5, 4, 6, 5))
#' res <- lsbclust(data = dat[[1]]$data, nclust = c(5, 4, 6, 5))
#' cfsim(res, dat[[1]])
#' 
#' ## Rows fixed, balanced classes
#' set.seed(2)
#' dat <- rlsbclust(ndata = 1, nobs = 100, size = c(10, 8), nclust = c(5, 4, 6, 5), 
#'                  fixed = "rows")
#' res <- lsbclust(data = dat[[1]]$data, nclust = c(5, 4, 6, 5), fixed = "rows")
#' cfsim(res, dat[[1]])
#' 
#' ## Rows fixed, unbalanced classes
#' set.seed(3)
#' dat <- rlsbclust(ndata = 1, nobs = 100, size = c(10, 8), nclust = c(5, 4, 6, 5), 
#'                  fixed = "columns", 
#'                  clustsize = list(NULL, NULL, c(40, 25, 15, 10, 5, 5), c(40, 25, 15, 10, 10)))
#' res <- lsbclust(data = dat[[1]]$data, nclust = c(5, 4, 6, 5), fixed = "columns")
#' cfsim(res, dat[[1]])
rlsbclust <- function(ndata = 50L, nobs, size, nclust, clustsize = NULL, delta = rep(1L, 4L), 
                      ndim = 2L, alpha = 0.5, fixed = c("none", "rows", "columns"), 
                      err_sd = 1, svmins = 1, 
                      svmax = 6) {
  
  cll <- match.call()
  
  ## Checks
  if (length(delta) != 4L)
    stop("Argument 'delta' should be of length 4.")
  if (any(!(delta %in% c(0L, 1L))))
    stop("Argument 'delta' should contain only 0's and 1's.")
  
  ## If nclust is of length one, expand it to length 4
  if (length(nclust) == 1) nclust <- rep(nclust, 4)
  if (length(nclust) != 4) stop("Argument 'nclust' should be either of length 1 or 4.")
  
  ## Expand svmins if needed
  if (length(svmins) == 1)
    svmins <- rep(svmins, ndim)
  
  ## Ensure that clustsize is length 4 (even if NULL)
  if (is.null(clustsize)) {
    clustsize <- as.list(floor(nobs / nclust))
    for (i in seq_along(clustsize)) {
      clustsize[[i]] <- rep(clustsize[[i]], nclust[i] - 1)
      clustsize[[i]] <- c(nobs - sum(clustsize[[i]]), clustsize[[i]])
    }
  }
  if (is.list(clustsize)) {
    ## Fill in NULL entries
    indNULL <- which(sapply(clustsize, is.null))
    if (length(indNULL)) {
      for (i in indNULL) {
        clustsize[[i]] <- floor(nobs / nclust[i])
        clustsize[[i]] <- rep(clustsize[[i]], nclust[i] - 1)
        clustsize[[i]] <- c(nobs - sum(clustsize[[i]]), clustsize[[i]])
      }
    }
    
    ## Do checks
    if (length(clustsize) != 4L) 
      stop("Argument 'clustsize' must be of length 4.")
    if (!isTRUE(all.equal(sapply(clustsize, sum), rep(nobs, length(clustsize))))) {
      cat("\n\nCluster sizes summed:", sapply(clustsize, sum), "\n\n")
      stop("The elements of 'clustsize' must sum to 'nobs'.")
    }
    if (any(sapply(clustsize, length) != nclust))
      stop("Number of elements in 'clustsize' does not correspond to 'nclust'.")
    if (any(sapply(clustsize, min) == 0))
      stop("Argument 'clustsize' should contain only positive integers.")
  }
 
  ## Process fixed
  fixed <- match.arg(tolower(fixed), choices = c("none", "rows", "columns"))
  
  ## Simulate cluster memberships with required probabilities (always generate all)
  cluster <- mapply(function(x, y) rep(seq_len(x), times = y)[sample.int(nobs)],
                    as.list(nclust), clustsize)
  colnames(cluster) <- c("overall", "rows", "columns", "interactions")
  
  ## Get cluster sizes
  osize <- table(cluster[, 1L])
  rsize <- table(cluster[, 2L])
  csize <- table(cluster[, 3L])
  isize <- table(cluster[, 4L])
  
  ## Simulate overall means, if necessary
  odelta <- delta[1] * delta[3] + delta[2] * delta[4] - delta[1] * delta[2]
  if (odelta) 
    omeans <- rnorm(nclust[1L])
  
  ## Simulate row means, if necessary, and centre
  if (delta[2L]) {
    rmeans <- replicate(n = nclust[2L], rnorm(size[1L]))
    if (delta[4L])
      rmeans <- scale(rmeans, scale = FALSE)
  }
  
  ## Simulate column means, if necessary
  if (delta[1L]) {
    cmeans <- replicate(n = nclust[3L], rnorm(size[2L]))  
    if (delta[3L])
      cmeans <- scale(cmeans, scale = FALSE)
  }
  
  ## Simulate interaction matrices
  ## Simulate orthornormal matrices
  if (fixed == "rows")
    Clst <- list(rorth(nrow = size[1L], ncol = ndim))
  else Clst <- replicate(nclust[4L], rorth(nrow = size[1L], ncol = ndim), 
                         simplify = FALSE)
  
  if (fixed == "columns")
    Dlst <- list(rorth(nrow = size[2L], ncol = ndim))
  else Dlst <- replicate(nclust[4L], rorth(nrow = size[2L], ncol = ndim), 
                         simplify = FALSE)
  
  ## Simulate singular values
  if (fixed == "none") {
    svlst <- simsv(nclust = nclust[4L], mins = svmins, max = svmax, ndim = ndim)
    # svlst <- replicate(nclust[4L], rev(cumsum(rlnorm(ndim, sdlog = sdlog, meanlog = meanlog))),
    #                    simplify = FALSE)
  } else {
    svlst <- simsv(nclust = 1L, mins = svmins, max = svmax, ndim = ndim)
    # svlst <- list(rev(cumsum(rlnorm(ndim, sdlog = sdlog, meanlog = meanlog))))
  }
  
  ## Multiply C and D with these
  Clst <- Map(function(x, y) x %*% diag(y ^ alpha), Clst, svlst)
  Dlst <- Map(function(x, y) x %*% diag(y ^ (1 - alpha)), Dlst, svlst)
  
  ## Centre these matrices, if required
  if (delta[1L]) {
    Dlst <- lapply(Dlst, function(x) scale(x, scale = FALSE))
  }
  if (delta[2L]) {
    Clst <- lapply(Clst, function(x) scale(x, scale = FALSE))
  }
  
  ## Function to construct interactions
  interfun.none <- function(i) tcrossprod(Clst[[i]], Dlst[[i]])
  interfun.rows <- function(i) tcrossprod(Clst[[1L]], Dlst[[i]])
  interfun.columns <- function(i) tcrossprod(Clst[[i]], Dlst[[1L]])
  interfun <- switch(fixed, none = interfun.none, 
                     rows = interfun.rows, 
                     columns = interfun.columns)
  
  ## Combine these into an array
  arr <- array(NA, dim = c(size, nobs))
  for (i in seq_len(nobs)) {
    arr[, , i] <- interfun(cluster[i, 4L])
    if (odelta)
      arr[, , i] <- arr[, , i] + omeans[cluster[i , 1L]]
    if (delta[2L])
      arr[, , i] <- arr[, , i] + rmeans[, cluster[i, 2L]] %o% rep(1L, size[2L]) 
    if (delta[1L])
      arr[, , i] <- arr[, , i] + rep(1L, size[1L]) %o% cmeans[, cluster[i, 3L]]
  }
  
  ## Add dimnames
  dimnames(arr) <- list(paste0("Row", seq_len(size[1L])), 
                        paste0("Col", seq_len(size[2L])), 
                        NULL)
  
  ## Create classes for different components (for \pkg{clue} compatibility)
  if (odelta) {
    omeans <- matrix(omeans, ncol = 1L)
    rownames(omeans) <- seq_len(nclust[1L])
    overall <- list(cluster = cluster[, "overall"], centers = omeans)
    class(overall) <- "lsbclust_sim_part"
  } else overall <- NULL
  if (delta[2L]) {
    rmeans <- t(rmeans)
    rownames(rmeans) <- seq_len(nclust[2L])
    rows <- list(cluster = cluster[, "rows"], centers = rmeans)
    class(rows) <- "lsbclust_sim_part"
  } else rows <- NULL
  if (delta[1L]) {
    cmeans <- t(cmeans)
    rownames(cmeans) <- seq_len(nclust[3L])
    columns <- list(cluster = cluster[, "columns"], centers = cmeans)
    class(columns) <- "lsbclust_sim_part"
  } else columns <- NULL
  interactions <- list(cluster = cluster[, "interactions"], C = Clst, D = Dlst,
                       svlst = svlst)
  class(interactions) <- "lsbclust_sim_part"
  
  ## Create the required number of data sets by adding random error
  simdata <- replicate(ndata, arr + array(rnorm(nobs * prod(size), sd = err_sd), 
                                          dim = c(size, nobs)), simplify = FALSE)
  
  ## convert to list of objects
  convfun <- function(x) {
    out <- list(data = x, actual = arr, cluster = cluster, overall = overall, 
                rows = rows, columns = columns, interactions = interactions, 
                sv = svlst, call = cll)
    class(out) <- c("lsbclust_sim", "list")  
    return(out)
  }
  
  out <- lapply(simdata, convfun)
  class(out) <- c("lsbclust_sim_collection", "list")  
  return(out)
}