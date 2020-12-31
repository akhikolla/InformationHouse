# fk_density computes the kernel estimated density from a univariate sample x.
# Arguments:
# x = sample points
# h = bandwidth. can be positive numeric, or one of 'Silverman' for Silverman's rule of thumb, or 'mlcv' for
#       maximum (pseudo) likelihood cross validation. default is Silverman's heuristic
# beta = vector of (positive) kernel coefficients. default is c(.25,.25)
# from,to = minimum,maximum of interval for grid/binned evaluation. default is min(x)-6h,max(x)+6h
# ngrid = number of grid points for evaluation. default is 1000
# nbin = number of bins if binning to be used for faster approximation. default is not to use binning (i.e. compute exact density estimates)
# x_eval = evaluation points. default is exact evaluation on a grid. if ngrid set to NULL then x_eval default is changed to the sample points themselves
#
# returns list with components
# $x = points at which evaluation of density was performed
# $y = estimated density values
# $h bandwidth used

fk_density <- function(x, h = 'Silverman', h_adjust = 1, beta = NULL, from = NULL, to = NULL, ngrid = 1000, nbin = NULL, x_eval = NULL){
  # check inputs
  if(!is.numeric(x)) stop('x must be a numeric vector')
  if(any(is.na(x))) stop('x cannot contain missing values')
  if(!is.null(beta) && (!is.vector(beta) || !is.numeric(beta))) stop('beta must be a numeric vector')
  if(!is.numeric(h_adjust) || length(h_adjust)>1 || h_adjust<0) stop('h_adjust must be a positive numeric')

  n <- length(x)

  # computation is slightly more stable for smaller magnitude sample, so centralise at zero
  mn <- mean(x)
  x <- x - mn
  if(!is.null(x_eval)) x_eval <- x_eval - mn
  if(!is.null(from)) from <- from - mn
  if(!is.null(to)) to <- to - mn

  # sample points must be sorted unless binning is used
  if(is.null(nbin) && is.unsorted(x)){
    o <- order(x)
    xo <-  x[o]
  }
  else xo <- x

  # set kernel coefficients, or normalise if provided
  if(is.null(beta)) beta <- c(.25, .25)
  else beta <- beta / norm_const_K(beta)

  # determine bandwidth value if actual value not provided
  if(!is.numeric(h)){
    if(h=='Silverman') h <- (8 * sqrt(pi) / 3 * roughness_K(beta) / var_K(beta)^2 / length(x))^.2 * sd(x)
    else if(h=='mlcv'){

      # mlcv bandwidth handled differently if binning is used
      if(!is.null(nbin)){
        frm <- min(x)
        tot <- max(x)
        wts <- sm_bin_wts(x, rep(1, n), nbin, frm, tot)
        xs <- seq(frm, tot, length = nbin)
        loo_cv <- function(h){
          hf <- ksum(xs, wts / (n - 1) / h, xs, h, beta, 1:nbin) - beta[1] / (n - 1) / h
          hf[hf < 1e-20] <- 1e-20 # leave-one-out values may be evaluated numerically to be zero, so buffer at small value
          sum(log(hf) * wts)
        }
        h <- suppressWarnings(optimise(loo_cv, sd(x) / n^.2 * c(.05, 5), maximum = TRUE)$maximum)
      }
      else{
        loo_cv <- function(h){
          hf <- ksum(xo, rep(1 / (n - 1) / h, n), xo, h, beta, 1:n) - beta[1] / (n - 1) / h
          hf[hf < 1e-20] <- 1e-20
          sum(log(hf))
        }
        h <- suppressWarnings(optimise(loo_cv, sd(x) / n^.2 * c(.05, 5), maximum = TRUE)$maximum)
      }
    }
    else stop('"h" must be positive numeric or one of "Silverman" or "mlcv"')
  }
  h <- h * h_adjust

  # evaluate density
  # three possibilities are binned, grid evaluation or neither
  if(!is.null(nbin)){ # binning

    if(!is.numeric(nbin)) stop('nbin must be positive natural number. Suggest at least 200')

    # establish range for grid of bin locations
    if(is.null(from)) from <- min(x) - 6 * h
    if(is.null(to)) to <- max(x) + 6 * h

    # compute bin values
    wts <- sm_bin_wts(x, rep(1, n), nbin, from, to)

    # set evaluation grid
    xs <- seq(from, to, length = nbin)

    # evaluation of density is handled differently if evaluation at the bin locations or specified set
    # of evaluation points
    if(is.null(x_eval)){

      # evaluate the density at bin locations
      y <- ksum(xs, wts, xs, h, beta, 1:nbin)

      # ensure normalisation is correct despite binning
      y <- y / sum(y) / (xs[2] - xs[1])
    }
    else{

      if(!is.numeric(x_eval)) stop('x_eval must be a numeric vector')

      # determine location of evaluation points in bins
      alloc <- cbin_alloc(x_eval, nbin, from, to)

      # evaluate density at evaluation points
      y <- ksum(xs, wts, x_eval, h, beta, alloc) / n / h
      xs <- x_eval
    }
  }
  else if(!is.null(ngrid) && is.null(x_eval)){ # exact evaluation on a grid

    if(!is.numeric(ngrid)) stop('ngrid must be a positive natural number')

    # establish range for grid
    if(is.null(from)) from <- xo[1] - 6 * h
    if(is.null(to)) to <- xo[n] + 6 * h

    # set evaluation points to grid locations
    xs <- seq(from, to, length = ngrid)

    # evaluate density with sorted sample at grid locations
    y <- ksum(xo, rep(1, n), xs, h, beta) / n / h
  }
  else{ # exact evaluation at non-grid locations
    if(is.null(x_eval)){ # default is evaluation at the sample points
      xs <- x
      xe <- xo
    }
    else{ # otherwise ensure evaluation points are sorted for fast kernel summing
      if(!is.numeric(x_eval)) stop('x_eval must be a numeric vector')
      xs <- x_eval
      if(is.unsorted(x_eval)) xe <- sort(x_eval)
      else xe <- x_eval
    }

    # evaluate density at evaluation points
    y <- ksum(xo, rep(1, n), xe, h, beta)[match(xs, xe)] / n / h
  }

  # in cases where x_eval contains values very far from actual sample, then numerical issues can make the outputs from c++ functions include NA or +- Inf. We set these to very small positive value

  y[is.na(y) | !is.finite(y)] <- 1e-300

  # return evaluation points (ensure to add the mean which was subtracted initially) and corresponding density values, plus bandwidth
  structure(list(x = xs + mn, y = y, h = h, call = match.call(), n = n), class = "fk_density")
}


# Methods for class fk_density

plot.fk_density <- function(x, main = NULL, ...){
  if(is.null(main)){
    main <- paste('fk_density(', names(x$call)[2], ' = ', as.character(x$call)[2], sep = '')
    if(length(x$call)>2) for(argnum in 3:length(x$call)) main <- paste(main, ', ', names(x$call)[argnum], ' = ', as.character(x$call)[argnum], sep = '')
    main <- paste(main, ')', sep = '')
  }
  plot(x$x, x$y, main = main, ylab = 'Density', xlab = paste('n = ', x$n, ', bandwidth = ', round(x$h, 4)), type = 'l', ...)
}

print.fk_density <- function(x, ...){
  cat('Call: \n \n')
  print(x$call)
  cat('\n')
  cat(paste('Data: x (', x$n, ' obs.); bandwidth = ', round(x$h, 4), ' \n \n', sep = ''))
}
