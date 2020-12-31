myoptim <- function(par, fn, gr,method,
                          lower, upper,
                          control, hessian,...)
{
  myfn <- function(par) {
    temp <- fn(par)
    if (!is.finite(temp)) temp <- 1.0e100
    return(temp)
  }
  if (is.null(control$trace)) verbose <- FALSE
  else verbose <- control$trace

    if (verbose) cat("using nlminb","\n")
  if (length(lower) < length(par)) lower <- rep(lower,length(par))
  if (length(upper) < length(par)) upper <- rep(upper,length(par))
  nlmcontrol <- control
  thenlm <- suppressWarnings(nlminb(par, fn,lower = lower, upper = upper, control=nlmcontrol))
   if (thenlm$convergence==0) theoptim <- list(par=thenlm$par,value=thenlm$objective,convergence=thenlm$convergence,message=thenlm$message)
  else {
    if (verbose) cat("using optim","\n")
    optimcontrol <- list(trace=ifelse(verbose,1,0),maxit=control$iter.max)
    theoptim <- suppressWarnings(optim(thenlm$par, myfn, method="L-BFGS-B",lower=lower,upper=upper,hessian=FALSE,control=optimcontrol))
  }
  if (theoptim$convergence!=0)  {
    if (verbose) cat("using constrOptim","\n")
    ui <- diag(length(par))[is.finite(lower),]
    ci <- lower[is.finite(lower)]
    theoptim <- suppressWarnings(constrOptim(par, fn, NULL , ui, ci, control=optimcontrol))
  }
  return(theoptim)
}