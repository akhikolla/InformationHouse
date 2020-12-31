# fk_sum returns kernel weighted sum, and/or kernel derivative weighted sum of terms omega based on locations x
# Arguments:
# x = vector sample points / kernel locations
# omega = vector of terms to be summed
# h = positive numeric bandwidth for kernel weights
# x_eval = points at which to evaluate sums. Default is evaluation at the sample points themselves
# beta = vector of kernel coefficients. Default is c(.25, .25)
# nbin = number of bins if binning approximation desired. Default is exact evaluation (nbin = NULL)
# type = one of 'ksum' for kernel weighted sum, 'dksum' for kernel derivative weighted sum,
#         or 'both' for both. Default is the kernel weighted sums (i.e. type is 'ksum')
#
# returns a vector of kernel or kernel derivative sums if type is 'ksum' or 'dksum', or a matrix
# with first column the kernel sums and second column the kernel derivative sums if type  is 'both'

fk_sum <- function(x, omega, h, x_eval = NULL, beta = c(.25, .25), nbin = NULL, type = "ksum"){
  # check inputs
  if(!is.numeric(x) || !is.numeric(omega) || length(x)!=length(omega)) stop('x and omega must be numeric vectors of the same length')
  if(any(is.na(c(x, omega)))) stop('x and omega cannot contain missing values')
  if(!is.vector(beta) || !is.numeric(beta)) stop('beta must be a numeric vector')
  if(!is.numeric(h) || h<0 || length(h)!=1) stop('h must be a positive numeric')

  n <- length(x)

  # computation is slightly more stable for smaller magnitude kernel locations, so centralise x at zero
  mn <- mean(x)
  x <- x - mn

  if(is.null(nbin)){

    # binning approximation not used, so sample needs to be sorted according to x values
    o <- order(x)
    xo <- x[o]
    omo <- omega[o]
    if(is.null(x_eval)){
      # evaluation points not supplied, so evaluate at x
      if(type == "ksum") ksum(xo, omo, x, h, beta, match(x, xo))
      else if(type=="dksum") dksum(xo, omo, x, h, beta, match(x, xo))
      else if(type=="both") kndksum(xo, omo, x, h, beta, match(x, xo))
    }
    else{
      if(!is.numeric(x_eval)) stop('x_eval must be a numeric vector')
      # evaluation points supplied, so subtract mean of x from these as well and sort, before computing sums
      x_eval <- x_eval - mn
      o_eval <- order(x_eval)
      x_eo <- x_eval[o_eval]
      if(type == "ksum") ksum(xo, omo, x_eo, h, beta)[match(x_eval, x_eo)]
      else if(type == "dksum") dksum(xo, omo, x_eo, h, beta)[match(x_eval, x_eo)]
      else if(type == "both") kndksum(xo, omo,x_eo, h, beta)[match(x_eval, x_eo),]
    }
  }
  else{ # binning approximation used so compute bin locations and binned omega values
    xo <- seq(min(x) - 1e-5, max(x) + 1e-5, length = nbin)
    omo <- sm_bin_wts(x, omega, nbin, xo[1], xo[nbin])
    if(is.null(x_eval)){
      # evaluation points not supplied, so evaluate at x
      if(type == "ksum") ksum(xo, omo, x, h, beta, cbin_alloc(x, nbin, xo[1], xo[nbin]))
      else if(type == "dksum") dksum(xo, omo, x, h, beta, cbin_alloc(x, nbin, xo[1], xo[nbin]))
      else if(type == "both") kndksum(xo, omo, x, h, beta, cbin_alloc(x, nbin, xo[1], xo[nbin]))
    }
    else{
      if(!is.numeric(x_eval)) stop('x_eval must be a numeric vector')
      # evaluation points supplied, so subtract mean of x from these as well, before computing sums
      x_eval <- x_eval - mn
      if(type == "ksum") ksum(xo, omo, x_eval, h, beta, cbin_alloc(x_eval, nbin, xo[1], xo[nbin]))
      else if(type == "dksum") dksum(xo, omo, x_eval, h, beta, cbin_alloc(x_eval, nbin, xo[1], xo[nbin]))
      else if(type == "both") kndksum(xo, omo, x_eval, h, beta, cbin_alloc(x_eval, nbin, xo[1], xo[nbin]))
    }
  }
}


# norm_const_K computes the normalising constant for a given vector of kernel coefficients.
# Arguments:
# beta = vector of kernel coefficients.
#
# returns a numeric C for which the kernel with coefficients beta/C has unit integral

norm_const_K <- function(beta){
  2 * sum(sapply(1:length(beta), function(k) beta[k] * factorial(k - 1)))
}

# roughness_K computes the squared L2 norm of the (normalised) kernel with coefficients beta. This is also
# sometimes referred to as the roughness of the kernel.
# Arguments:
# beta = vector of kernel coefficients.

roughness_K <- function(beta){
  beta_use <- beta / norm_const_K(beta)
  betakj <- beta_use %*% t(beta_use)
  kpj <- log((2^(0:(length(beta) - 1))) %*% t(2^(0:(length(beta) - 1))), base = 2)
  sum(betakj / 2^kpj * factorial(kpj))
}

# var_K computes the variance of the random variable with kernel coefficients beta / norm_const_K(beta).
# Arguments:
# beta = vector of kernel coefficients.

var_K <- function(beta){
  beta_use <- beta / norm_const_K(beta)
  2 * sum(sapply(1:length(beta), function(k) beta_use[k] * factorial(k + 1)))
}

# h_Gauss_to_K transforms a bandwidth for use with the Gaussian kernel to one which is appropriate for use with
# the kernel with coefficients beta.
# Arguments:
# h = positive numeric bandwidth for Gaussian kernel
# beta = vector of kernel coefficients.

h_Gauss_to_K <- function(h, beta){
  h * (roughness_K(beta) / var_K(beta)^2 * 8 * sqrt(pi) / 4)^(1 / 5)
}

# plot_kernel produces a visualisation of the shape of a kernel with coefficients beta.
# Arguments:
# beta = vector of kernel coefficients
# type = (optional) character. The plot type passed to the base plot() function
# ... = (optional) additional arguments to be passed to the base plot() function

plot_kernel <- function(beta, type = 'l', ...){
  sig <- sqrt(var_K(beta))
  nc <- norm_const_K(beta)
  xs <- seq(-4, 4, length = 500)
  ys <- rowSums(sweep(sapply(0:(length(beta) - 1), function(k) abs(xs * sig)^k), 2, beta, '*') * exp(-abs(xs * sig)))
  plot(xs, ys / nc * sig, xlab = 'x', ylab = 'K(x)', type = type, ...)
}

# norm_K computes the L2 norm of the (normalised) kernel with coefficients beta.
# Arguments:
# beta = vector of kernel coefficients.

norm_K <- function(beta) sqrt(roughness_K(beta))

# h_K_to_Gauss transforms a bandwidth for use with the kernel with coefficients beta to one which is appropriate for use
# with the Gaussian kernel.
# Arguments:
# h = positive numeric bandwidth for Gaussian kernel
# beta = vector of kernel coefficients.

h_K_to_Gauss <- function(h, beta){
  h / (roughness_K(beta) / var_K(beta)^2 * 8 * sqrt(pi) / 4)^(1 / 5)
}
