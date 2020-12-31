# f_ica evaluates the objective function for ICA based on (negative) sample entropy estimation
# Arguments:
# w = projection vector. That is, function returns hat_entropy(Xw)
# X = data matrix
# h = bandwidth for density estimation
# betas = vector of kernel coefficients
# nbin = number of bins if binning approximation to be used. default is exact evaluation
#
# returns numeric, the value of the objective

f_ica <- function(w, X, h, betas, nbin = NULL){

  # compute projected data onto normalised w
  x <- c(X %*% w / sqrt(sum(w^2)))
  n <- length(x)


  if(is.null(nbin)){ # no binning, so exact density estimate computed at sample points
    # fast kernels sums rely on ordered sample
    o <- order(x)
    xo <- x[o]

    # compute (scaled) density values at the sample and buffer zeroes
    hf <- ksum(xo, rep(1, n), xo, h, betas, 1:n)
    hf[which(hf < 1e-300)] <- 1e-300

    # return negative of sample entropy
    sum(log(hf))
  }
  else{ # binning approximation used
    # set up grid of bin locations
    m <- min(x) - 1e-5
    M <- max(x) + 1e-5
    xs <- seq(m, M, length = nbin)

    # compute bin counts
    ynew <- bin_wts(x, rep(1, n), nbin, m, M)

    # approximate density at bin locations and buffer zeroes
    hf <- ksum(xs, ynew, xs, h, betas, 1:nbin)
    hf[which(hf < 1e-300)] <- 1e-300

    # return approximate negative sample entropy
    sum(log(hf) * ynew)
  }
}


# df_ica computes the gradient of the ICA objective based on (negative) sample entropy
# of projected data
# Arguments:
# w = projection vector. That is, function returns d (hat_entropy(Xw))/ d w
# X = data matrix
# h = bandwidth for density estimation
# betas = vector of kernel coefficients
# nbin = number of bins if binning approximation to be used. default is exact evaluation
#
# returns numeric vector, the gradient of the objective

df_ica <- function(w, X, h, betas, nbin = NULL){

  # compute projected data onto normalised projection vector
  nw <- sqrt(sum(w^2))
  x <- c(X %*% w / nw)
  n <- length(x)
  if(is.null(nbin)){ # no binning, so exact density estimate computed at sample points

    # fast kernel sums rely on ordered sample
    o <- order(x)
    xo <- x[o]

    # compute (scaled) density and buffer zeroes
    hf <- ksum(xo, rep(1, n), xo, h, betas, 1:n)
    hf[which(hf < 1e-300)] <- 1e-300

    # compute other components of partial derivatives of -entropy w.r.t. projected points
    dhf <- -dksum(xo, rep(1, n), xo, h, betas, 1:n) / h
    d2 <- -dksum(xo, 1 / hf, xo, h, betas, 1:n) / h

    # dp represents the partial derivatives of -entropy w.r.t. projected points
    dp <- (d2 + dhf / hf)

    # chain rule product with dp/dw gives total gradient of the objective w.r.t. w
    dp %*% (X[o,] / nw - xo %*% t(w) / nw^2)
  }
  else{ # binned approximation used

    # setup bin locations and bin weights
    m <- min(x) - 1e-5
    M <- max(x) + 1e-5
    ynew <- bin_wts(x, rep(1, n), nbin, m, M)
    xs <- seq(m, M, length = nbin)

    # determine allocation of projected points to bins
    alloc <- cbin_alloc(x, nbin, m, M)

    # compute partial derivatives w.r.t. bin locations, similar to above
    hf <- ksum(xs, ynew, xs, h, betas, 1:nbin)
    hf[which(hf < 1e-300)] <- 1e-300
    dhf <- - dksum(xs, ynew, xs, h, betas, 1:nbin) / h
    ynew <- bin_wts(xs, 1 / hf * ynew, nbin, m, M)
    d2 <- - dksum(xs, ynew, xs, h, betas, 1:nbin) / h
    dp <- (d2 + dhf / hf)

    # compute gradient of objective using chain rule product, as before.
    # approximate partial derivative w.r.t. projected point x[i] is dp[alloc[i]]
    dp[alloc] %*% (X / nw - x %*% t(w) / nw^2)
  }
}


# fk_ICA performs independent component analysis based on kernel based estimate of
# sample entropy.
# Arguments:
# X = data matrix from which to extract components
# ncomp = number of independent components. Currently ICs are computed from first ncomp
#         principal components, as in fastICA. For solution from more flexible optimisation
#         use fk_ICA(X, ncomp + k) for some k > 0 and take only first ncomp terms from output.
# beta = kernel coefficients. default is c(.25, .25)
# hmult = bandwidth multiplier. bandwidth uses is Silverman's rule of thumb multiplied
#         by hmult.
# it = number of iterations in optimisation.
# nbin = number of bins if binning approximation is to be used. Default is based on exact
#         density estimate.
#
# returns list with components
# $X = data matrix passed to function
# $K = whitening matrix
# $W = unmixing matrix for whitened data
# $S = independent components, i.e., S = sweep(X, 2, colMeans(X), '-')%*%K%*%W


fk_ICA <- function(X, ncomp = 1, beta = c(.25, .25), hmult = 1.5, it = 20, nbin = NULL){
  # check inputs
  if(!is.matrix(X)) stop('X must be a numeric matrix')
  if(any(is.na(X))) stop('X cannot contain missing values')
  if(!is.numeric(ncomp)) stop('ncomp must be a positive natural number')
  if(!is.vector(beta) || !is.numeric(beta)) stop('beta must be a numeric vector')
  if(!is.numeric(hmult) || length(hmult)>1 || hmult<0) stop('hmult must be a positive numeric')
  if(!is.numeric(it)) stop('the number of iterations (it) must be a positive natural number')
  if(!is.null(nbin) && !is.numeric(nbin)) stop('nbin must be a positive natural number')

  dim <- ncomp
  n <- nrow(X)

  # Independent components computed from whitened data
  W <- whiten(X, dim)
  Y <- W$Y

  # V stores independent component projections based on whitened data
  V <- matrix(0, dim, dim)

  # Determine ICs iteratively. Only ncomp-1 need to be found, since the final
  # IC will be any vector in the null space of the components so far found
  for(i in 1:(dim-1)){

    # v0 is the projection vector to be optimised
    # all those beyond the first must be made orthogonal to
    # those found before
    v0 <- numeric(dim)
    v0[i] <- 1
    if(i > 1){
      for(j in 1:(i - 1)){
        v0 <- v0 - (v0 %*% V[,j])[1] * V[,j]
        v0 <- v0 / sqrt(sum(v0^2))
      }
    }

    # compute bandwidth using hmult*h(Silverman)
    h <- hmult * (norm_K(beta)^2 / var_K(beta)^2 * 8 * sqrt(pi) / 3 / nrow(X))^(1 / 5)

    # the following is a very simple gradient ascent. optim could also be used, but
    # this has been found to work efficiently and effectively for ICA
    for(iter in 1:it){
      fval <- f_ica(v0, Y, h, beta, nbin) # function value at current iterate
      dir <- c(df_ica(v0, Y, h, beta, nbin)) # search direction
      ndir <- sqrt(sum(dir^2))
      dir <- dir / ndir
      stp <- .5 # step size for backtracking line-search
      repeat{ # line-search
        fnew <- f_ica(v0 + stp * dir, Y, h, beta, nbin)
        if(fnew > (fval / .99999)) break
        else stp <- stp * .5
        if(stp < 1e-9) break
      }
      v0 <- v0 + stp * dir
      v0 <- v0 / sqrt(sum(v0^2))
      if(stp < 1e-9) break
    }
    v0 <- v0 / sqrt(sum(v0^2))

    # if this is not the first component, then ensure orthogonality
    # is retained
    if(i > 1){
      for(j in 1:(i - 1)){
        v0 <- v0 - (v0 %*% V[,j])[1] * V[,j]
        v0 <- v0 / sqrt(sum(v0^2))
      }
    }

    # project data into null-space of the component just found to
    # ensure it isn't rediscovered. This is a common approach in
    # projection pursuit to avoid having to use a constrained
    # optimisation implementation
    Y <- Y - Y %*% v0 %*% t(v0)
    V[,i] <- v0
  }
  V[,dim] <- Null(V[,1:(dim-1)])
  structure(list(X = X, K = W$E$vectors[,1:dim] %*% diag(1 / sqrt(W$E$values[1:dim])), W = V, S = W$Y %*% V, call = match.call()), class = 'fk_ICA')
}



# whiten performs data whitening (i.e., "removal" of scale and location)
# Argments:
# X = data matrix
# ncomp = number of whitened components desired
#
# returns list with components
# $E = whitening matrix
# $Y = whitened data

whiten <- function(X, ncomp){

  # if data are relatively high dimensional then use iterative eigen-solver
  # from rARPACK. Otherwise use standard function eigen()
  if(ncol(X) > 300) E <- eigs_sym(cov(X), ncomp)
  else E <- eigen(cov(X))

  # if not sufficient rank, use pseudo-inverse components
  E$values[which(E$values < 1e-10)] <- Inf
  Y <- t(t(X) - colMeans(X)) %*% E$vectors[, 1:ncomp] %*% diag(1 / sqrt(E$values[1:ncomp]))
  list(E = E, Y = Y)
}

plot.fk_ICA <- function(x, ...){
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  for(comp in 1:ncol(x$S)){
    readline("Press [Enter] to view the next component plot. Press [Esc] to exit.")
    plot(x$S[,comp], main = paste('Plot of component ', comp, sep = ''), xlab = 'index', ylab = 'Component value', ...)
    plot(fk_density(x$S[,comp]), main = paste('Esimated density of component ', comp, sep = ''), ...)
  }
  par(op)
}

print.fk_ICA <- function(x, ...){
  cat('Call: \n \n')
  print(x$call)
  cat('\n')
  cat(paste("Data: X (", nrow(x$X), " obs. in ", ncol(x$X), " dimensions); \n", sep = ''))
  cat(paste("Number of estimated source components = ", ncol(x$S),  "\n", sep = ''))
}
