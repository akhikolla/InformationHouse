fun_est_parallel <- function(max.theta, tilde_b, arglist, cl) {

  raw.max.theta <- max.theta1 <- 10^ceiling(log(max.theta, base = 10))
  theta.dec <- 1
  seq.theta <- ceiling(seq(from = max.theta1, to = 0, by = - max.theta1/10))
  stop.rule <- 1
  depth <- 1

  while (stop.rule >= 1) {
    several.fits <- fun_penSurvIC(theta = seq.theta[depth], tilde_b, arglist, tol = 0.01)$b
    if(any(several.fits != 0)) {
      if(depth == 1) stop("Please increase 'max.theta1'.")
      max.theta1 <- seq.theta[depth - 1]
      min.theta <- seq.theta[depth]
      seq.theta <- seq(from = max.theta1, to = min.theta, by = -raw.max.theta/(10^(theta.dec+1)))
      stop.rule <- raw.max.theta/(10^(theta.dec + 1))
      theta.dec <- theta.dec + 1
      depth <- 1
    } else {
      depth <- depth + 1
    }
  }

  upper.theta <- max.theta1
  e <- 0.0001;  rr <- 100
  lower.theta <- e * upper.theta
  set.theta <- upper.theta*(lower.theta/upper.theta)^((0:rr)/rr)
  set.theta <- matrix(set.theta, ncol = 1)

  if (is.null(cl)) {
    bic_b <- t(apply(set.theta, 1, fun_bic, tilde_b = tilde_b, arglist = arglist, tol = 0.01))
  } else {
    if (inherits(cl, "cluster")) {
      parallel_fun <- if (isTRUE(getOption("pboptions")$use_lb)) parLapplyLB else parLapply
      bic_b0 <- parallel_fun(cl, set.theta, fun_bic, tilde_b = tilde_b, arglist = arglist, tol = 0.01)
    } else {
      bic_b0 <- mclapply(set.theta, fun_bic, tilde_b = tilde_b, arglist = arglist, tol = 0.01)
    }
    bic_b <- t(sapply(bic_b0, function(x) x))
  }


  o.bic <- which.min(bic_b[,1])
  final.bic <- bic_b[o.bic, 1]
  final.theta <- set.theta[o.bic]
  results <- fun_penSurvIC(final.theta, tilde_b, arglist)
  results$bic <- final.bic

  return(results)

}
