#' T3Clusf: Tucker3 Fuzzy Cluster Analysis
#' 
#' This is an implementation of the T3Clusf algorithm of Rocci & Vichi (2005).
#' 
#' @param X Three-way data array, with no missing values.
#' @param Q Integer giving the number of dimensions required for mode B (variables).
#' This is the first mode of the array, excluding the mode clustered over (see \code{margin}).
#' @param R Integer giving the number of dimensions required for mode C (occasions). 
#' This is the second mode of the array, excluding the mode clustered over (see \code{margin}).
#' @param G Integer giving the number of clusters required.
#' @param margin Integer giving the margin of the array to cluster over. The remaining two
#' modes, in the original order, corresponds to \code{Q} and \code{R}.
#' @param alpha Numeric value giving the fuzziness parameter.
#' @param eps Small numeric value giving the empirical convergence threshold.
#' @param maxit Integer giving the maximum number of iterations allowed.
#' @param verbose Integer giving the number of iterations after which the loss values are printed.
#' @param nstart Integer giving the number of random starts required.
#' @param parallel Logical indicating whether to parallelize over random starts if 
#' \code{nstart > 1}.
#' @param mc.cores Argument passed to \code{\link{makeCluster}}.
#' @param minsize Integer giving the minimum size of cluster to uphold when reinitializing
#' empty clusters.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' @references 
#' Rocci, R., & Vichi, M. (2005). \emph{Three-mode component analysis with crisp or fuzzy partition of units}. 
#' Psychometrika, 70(4), 715-736.
#' @examples 
#' data("dcars")
#' set.seed(13)
#' res <- T3Clusf(X = carray(dcars), Q = 3, R = 2, G = 3, alpha = 1)
#' 
T3Clusf <- function(X, Q, R = Q, G = 2, margin = 3L, alpha = 1, eps = 1e-8, maxit = 100L,
                    verbose = 1, nstart = 1L, parallel = TRUE, 
                    mc.cores = detectCores() - 1L, minsize  = 3L) {
  
  ## Recurse if multiple starts needed
  if (nstart > 1L) {
    if (parallel) {
      if (.Platform$OS.type == "windows") {
        cl <- makeCluster(mc.cores, type = "PSOCK")
      } else {
        cl <- makeCluster(mc.cores, type = "FORK")
      }
      registerDoParallel(cl)
      out <- foreach (i = seq_len(nstart)) %dopar%  {
        T3Clusf(X = X, Q = Q, R = R, G = G, margin = margin, alpha = alpha, eps = eps, 
                maxit = maxit, verbose = verbose, nstart = 1L)
      }
      stopCluster(cl)
    } else {
      out <- replicate(nstart, T3Clusf(X = X, Q = Q, R = R, G = G, margin = margin, 
                                          alpha = alpha, eps = eps, maxit = maxit, 
                                          verbose = verbose, nstart = 1L), simplify = FALSE)
    }
    
    ## Determine the start with the best loss
    allloss <- sapply(out, "[[", "minloss")
    wmin <- which.min(allloss)
    
    ## Return only this start, with loss added
    out <- out[[wmin]]
    out$allloss <- allloss
    return(out)
  }
  
  ## Call
  cll <- match.call()
  
  ## Permute array so that clustering is over FIRST way
  if (!is(X, "array")) stop("'X' must be an array.")
  if (length(dim(X)) != 3) stop("'X' must be a three-dimensional array.")
  if (!(margin %in% seq_len(3)))  stop("'margin' must be either 1, 2, or 3.")
  perm <- c(margin, seq_len(3)[-margin])
  X <- aperm.default(X, perm = perm)
  dims <- dim(X)
  Xmat <- matrix(X, nrow = dims[1], ncol = prod(dims[-1]))
  
  ## Parameter checks
  if (any(is.na(X))) stop("'X' contains missing values.")
  if (alpha < 1) stop("'alpha' must be larger than or equal to 1.")
  
  ## Starting values
  ## Orthogonal B
  B <- matrix(rnorm(dims[2] * Q), nrow = dims[2], ncol = Q)
  B <- qr(B)
  B <- qr.Q(B)
  ## Orthogonal C
  C <- matrix(rnorm(dims[3] * R), nrow = dims[3], ncol = R)
  C <- qr(C)
  C <- qr.Q(C)
  ## Membership matrix U
  if (alpha == 1) {
    U <- sample.int(n = G, size = dims[1], replace = TRUE)
    U <- diag(G)[U, ]
  } else {
    U <- matrix(runif(dims[1] * G), nrow = dims[1], ncol = G)
    U <- diag(1 / rowSums(U)) %*% U
  }
  ## Starting values of Xmns (KJ x G) and Y (G x QR)
  Xmns <- apply(U^alpha, 2, function(z) apply(z * X, 2:3, sum) / sum(z))
  Y <- crossprod(Xmns, kronecker(C, B))
  
  ## Denominator for standardizing the loss
  maxloss <- sum(X^2)
  
  ## Monitor loss
  loss <- rep(NA, maxit)
  iter <- 0L
  
  ## Monitor empty clusters
  nempty <- rep(0L, maxit)
  iterempty <- rep(FALSE, maxit)
  
  ## Iterate
  while (iter < maxit) {
    
    iter <- iter + 1L
    
    ## Update U
    CBY <- tcrossprod(kronecker(C, B), Y)
    dists <- apply(Xmat, 1, function(x, y) colSums((x - CBY)^2))
    if (alpha == 1) {
      U <- diag(G)[apply(dists, 2, which.min), ]
    } else {
      dists <- dists^(-1/(alpha - 1))
      U <- t(dists / matrix(colSums(dists), nrow = G, ncol = dims[1], byrow = TRUE))
    }
    
    ## Check for empty clusters if alpha == 1
    if (alpha == 1) {
      clustsize <- colSums(U)
      zeroclass <- seq_len(G)[clustsize == 0]
      if (length(zeroclass) > 0) {
        
        ## Monitor empty clusters
        iterempty[iter] <- TRUE
        nempty[iter] <- sum(clustsize == 0)
        
        ## Reinitialize empty cluster(s) with worst fitting observation(s)
        ## Cluster membership and distances
        clustmem <- apply(dists, 2L, which.min)
        clustdists <- dists[cbind(clustmem, seq_len(dims[1]))]
        
        ## Order from worst to best-fitting
        ord <- order(clustdists, decreasing = TRUE)
        
        ## Remove objects from classes smaller than or equal to minsize from consideration
        smallclasses <- which(clustsize <= minsize)
        ord <- setdiff(ord, which(clustmem %in% smallclasses))
        
        ## Move worst fitting observations in consideration set to empty cluster(s)
        id <- ord[seq_len(sum(clustsize == 0))]
        U[id, ] <- 0L
        U[cbind(id, zeroclass)] <- 1
        
        ## Print message
        message("T3Clusf: ", sum(clustsize == 0), " empty cluster(s) re-iniatilized.")
        
        ## Update clustsize
        clustsize <- colSums(U)
      }
    }
    
    ## Update Xmns and Y
    Xmns <- apply(U^alpha, 2, function(z) apply(z * X, 2:3, sum) / sum(z))
    Y <- crossprod(Xmns, kronecker(C, B))

    ## Update B
    omega <- colSums(U^alpha)
    CC <- tcrossprod(C)
    bmat <- matrix(0, nrow = dims[2], ncol = dims[2])
    for (i in seq_len(G)) {
      matX <- matrix(Xmns[, i], nrow = dims[2], ncol = dims[3])
      bmat <- bmat + omega[i] * tcrossprod(matX %*% CC, matX)
    }
    B <- eigen(bmat)$vectors[, seq_len(Q)]
    Y <- crossprod(Xmns, kronecker(C, B))
    
    ## Update C
    BB <- tcrossprod(B)
    cmat <- matrix(0, nrow = dims[3], ncol = dims[3])
    for (i in seq_len(G)) {
      matX <- matrix(Xmns[, i], nrow = dims[2], ncol = dims[3])
      cmat <- cmat + omega[i] * crossprod(matX, BB) %*% matX
    }
    C <- eigen(cmat)$vectors[, seq_len(R)]
    Y <- crossprod(Xmns, kronecker(C, B))
    
    ## Calculate the loss
    CC <- tcrossprod(C)
    CCBB <- kronecker(CC, BB)
    target <- CCBB %*% Xmns
    losscomps <- matrix(nrow = dims[1], ncol = G)
    for (i in seq_len(dims[1])) {
      losscomps[i, ] <- colSums((matrix(Xmat[i, ], nrow = prod(dims[-1]), ncol = G) - target)^2)
    }
    indloss <- rowSums(U^alpha * losscomps)
    loss[iter] <- sum(indloss) / maxloss
    
    ## Monitor convergence
    if (verbose > 0 && iter %% verbose == 0)
      cat(sprintf(paste0("%", nchar(as.character(maxit)), "d"), iter), 
          "| Loss =", sprintf("%5f", loss[iter]), "\n")
    
    ## Check convergence
    if (iter > 2) {
      if (loss[iter] > loss[iter - 1]) 
        warning("The loss increased in iteration ", iter)
      if (1 - loss[iter] / loss[iter - 1] < eps) break
    }
  }
  
  ## Add rownames
  rownames(B) <- dimnames(X)[[2]]
  rownames(C) <- dimnames(X)[[3]]
  
  ## Reshape Y into array
  Yarr <- array(t(Y), dim = c(Q, R, G))
  
  ## Calculate cluster means
  means <- alply(Yarr, 3L, function(x) B %*% tcrossprod(x, C))
  
  ## Create output list
  out <- list(G = U, B = B, C = C, H = Yarr, cluster = apply(U, 1, which.max), 
              means = means, iter = iter, loss = loss[seq_len(iter)], 
              minloss = loss[iter], sizes = colSums(U), 
              nempty = nempty, iterempty = iterempty, call = cll)
  
  ## Construct fitted data
  if (alpha == 1L) {
    fitted <- array(unlist(means[out$cluster]), dim = dim(X)[order(perm)])  
  } else {
    fitted <- array(NA, dim = dim(X))
    for (i in seq_len(dims[1L]))
      fitted[i, , ] <- Reduce("+", Map("*", U[i, ], means))
    fitted <- aperm(fitted, perm = order(perm))
  }
  dimnames(fitted) <- dimnames(X)[order(perm)]
  out$fitted <- fitted
  out$alpha <- alpha
  
  ## Make S3 class
  class(out) <- c("T3Clusf", "list")
  return(out)
}