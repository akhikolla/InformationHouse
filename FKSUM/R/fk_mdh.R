# fk_fmdh computes the objective function for minimum density projection pursuit
# That is, the (penalised) integrated density on the minimum density hyperplane
# orthogonal to the projection vector
# Arguments:
# w = projection vector
# X = data matrix. Must have mean zero!
# h = bandwidth used in density estimate
# beta = vector of kernel coefficients
# alpha = width of "feasible region". That is, penalty only applied outside (-al*sd, al*sd)
#         where sd is the standard deviation of the projected data
# C = strength of penalty. Higher values better approximate constraint, but very high values
#         can affect numerical stability
# n_eval = number of evaluation points in initial search for minimum density hyperplane
#
# returns numeric, the value of the objective function


fk_fmdh <- function(w, X, h, beta, alpha, C, n_eval = 100){

  # compute projected data onto normalised projection vector
  x <- X %*%w / sqrt(sum(w^2))
  s <- sd(x)
  n <- length(x)

  # compute value of the objective using Rcpp function fk_md
  fk_md(x, rep(1, n), 1.2 * seq(-alpha * s, alpha * s, length = n_eval), h, beta ,alpha , C) / n / h
}

# fk_dfmdh computes the gradient of the objective function for minimum density projection pursuit
# Arguments:
# w = projection vector
# X = data matrix. Must have mean zero!
# h = bandwidth used in density estimate
# beta = vector of kernel coefficients
# alpha = width of "feasible region". That is, penalty only applied outside (-al*sd, al*sd)
#         where sd is the standard deviation of the projected data
# C = strength of penalty. Higher values better approximate constraint, but very high values
#         can affect numerical stability
# n_eval = number of evaluation points in initial search for minimum density hyperplane
#
# returns numeric vector, the gradient of the objective function


fk_dfmdh <- function(w, X, h, beta, alpha, C, n_eval = 100){

  # compute projected data onto normalised projection vector
  nw <- sqrt(sum(w^2))
  x <- X %*% w / nw
  s <- sd(x)
  n <- length(x)

  # compute partial derivatives of objective w.r.t. projected points using Rcpp function fk_md_dp
  dp <- fk_md_dp(x, rep(1, n), 1.2 * seq(-alpha * s, alpha * s, length = n_eval), h, beta, alpha, C) / n / h^2

  # return chain rule product of dp %*% (derivative of projected points w.r.t. w)
  c(dp %*% (X / nw - x %*% t(w) / nw^2))
}



# fk_mdh computes the minimum density hyperplane for the purpose of (two-way) clustering of data X
# Arguments:
# X = data matrix.
# v0 = initial projection vector. default is the first principal component projection
# hmult = bandwidth multiplier. default bandwidth is hmult multiplied by Silverman's rule of thumb
#         computed from data projected on initial projection vector, v0
# beta = vector of kernel coefficients. default is c(.25, .25)
# alphamax = maximum width of "feasible region". Effectively enforces the minimum density hyperplane to
#         intersect the ellipsoid v'Sv < alphamax^2 in the space with origin equal to the mean of the data,
#         where S is the covariance of the data.
#
# returns list with components
# $v = optimal projection vector
# $b = location of optimal hyperplane along v

fk_mdh <- function(X, v0 = NULL, hmult = 1, beta = c(.25, .25), alphamax = 1){
  # check inputs
  if(!is.matrix(X)) stop('X must be a numeric matrix')
  if(any(is.na(X))) stop('X cannot contain missing values')
  if(!is.vector(beta) || !is.numeric(beta)) stop('beta must be a numeric vector')
  if(!is.numeric(hmult) || length(hmult)>1 || hmult<0) stop('hmult must be a positive numeric')
  if(!is.numeric(alphamax) || length(alphamax)>1 || alphamax<0) stop('alphamax must be a positive numeric')
  if(!is.null(v0) && (!is.vector(v0) || length(v0)!=ncol(X))) stop('v0 must be a numeric vector of length ncol(X)')

  n <- nrow(X)

  # centralise data since objective and gradient computations require zero mean data
  mn <- colMeans(X)
  X <- sweep(X, 2, mn, '-')

  # initialise projection vector
  if(is.null(v0)){
    if(ncol(X) > 300) v0 <- eigs_sym(cov(X), 1)$vectors[,1] # if high dimensional then use faster eigen-solver
    else v0 <- eigen(cov(X))$vectors[,1] # otherwise use standard eigen-solver
  }
  else v0 <- v0/sqrt(sum(v0^2))

  # compute bandwidth for density estimation
  h <- hmult * (8 * sqrt(pi) * roughness_K(beta) / var_K(beta)^2 / 3 / n)^.2 * sd(X %*% v0)

  # setup parameters of optimal solution since optimum may not be the final solution if
  # by increasing alpha the solution becomes "invalid" (doesn't separate modes of the density)
  alpha_opt <- 0
  vopt <- v0
  alpha <- 0
  b <- 0
  C <- 100000 # strenght of penalty term. optimisation is fairly insensitive to this value

  # increase size of feasible region until the maximum is reached, and return final valid solution found
  while(alpha < alphamax){

    # obtain optimal projection vector using BFGS algorithm and normalise
    v0 <- optim(v0, fk_fmdh, fk_dfmdh, X, h, beta, alpha, C, method = 'BFGS', control = list(maxit=50))$par
    v0 <- v0 / sqrt(sum(v0^2))

    # compute projected data onto optimal vector
    x <- X %*% v0
    s <- sd(x)

    # check if the solution is "valid", i.e., separates the modes of the projected density
    if(fk_is_minim_md(x, rep(1, n), 2 * seq(-alpha * s, alpha * s, length = 1000), h, beta, alpha, C)){
      # if the solution is valid then update the solution to be returned
      vopt <- v0
      alphaopt <- alpha
      b <- fk_md_b(x, rep(1, n), 1.1 * seq(-alpha * s, alpha * s, length = 1000), h, beta, alpha, C) + mn %*% v0
    }

    # increase size of feasible region and continue
    alpha <- alpha + .1
  }
  structure(list(v = v0, b = b, X = sweep(X, 2, mn, '+'), h = h, call = match.call()), class = 'fk_mdh')
}


# Methods for class fk_mdh

plot.fk_mdh <- function(x, ...){
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  x2 <- x$X - x$X%*%x$v%*%t(x$v)
  v2 <- eigen(cov(x2))$vectors[,1]
  plot(x$X%*%x$v, x$X%*%v2, main = 'Minimum density hyperplane', xlab = "Optimal projection", ylab = "First PC in Nullspace", ...)
  abline(v = x$b, col = 2, lwd = 2)
  plot(fk_density(x$X%*%x$v, h = x$h), main = 'Estimated density on optimal projection', ...)
  abline(v = x$b, col = 2, lwd = 2)
  par(op)
}

print.fk_mdh <- function(x, ...){
  cat('Call: \n \n')
  print(x$call)
  cat('\n')
  cat(paste("Data: X (", nrow(x$X), " obs. in ", ncol(x$X), " dimensions); \n", sep = ''))
}
