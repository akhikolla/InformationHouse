robmixglm <-
 function(formula,family=c("gaussian","binomial","poisson","gamma","truncpoisson","nbinom"),data,offset=NULL,quadpoints=21,
          notrials=50,EMTol=1.0e-4, cores = max(detectCores() - 1, 1), verbose=FALSE) {
      
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  call <- match.call()
  
  if (missing(family)) family <- "gaussian"
  
  if (!(family %in% c("gaussian","binomial","poisson","gamma","truncpoisson","nbinom")))
    stop("Valid families are gaussian, binomial, poisson, gamma, truncpoisson, nbinom.\n")
  
  if (missing(data)) data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","offset"), names(mf), 0L)
  mf <- mf[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  Y <- model.response(mf, "any")

  X <- model.matrix(mt,mf)

  offset <- model.extract(mf,"offset")

  if(is.null(offset)) offset <- rep(0.0,dim(X)[1])
    
  if (family=="binomial")  {
    if (!inherits(Y,"matrix")) stop("Binomial data must be in success failure form.")
    if (dim(Y)[2]!=2) stop("Binomial data must be in success failure form.")
    if (any(Y[,1] <0) | any(Y[,2] <0)) stop("Binomial data must be positive.")
    if (any(!is.wholenumber(Y[,1]) | !is.wholenumber(Y[,2]))) stop("Binomial data must be integers.")
   } else if ((class(Y)!="numeric") & (class(Y)!="integer")) stop("Data must be a single column") 
  
  if (family=="poisson") {
    if (any(Y <0)) stop("Poisson data must be positive.")
    if (any(!is.wholenumber(Y))) stop("Poisson data must be integers.")
  }
  
  if (family=="gamma") {
    if (any(Y <0)) stop("Gamma data must be positive.")
  }
  
  if (family=="truncpoisson") {
    if (any(Y<=0)) stop("Truncated poisson data must be positive non-zero.")
    if (any(!is.wholenumber(Y))) stop("Truncated poisson data must be integers.")
  }

  if (family=="nbinom") {
    if (any(Y <0)) stop("Negative binomial data must be positive.")
    if (any(!is.wholenumber(Y))) stop("Negative binomial data must be integers.")
  }
  
  if (verbose & (cores>1)){
    warning("Setting cores=1 to allow Verbose output")
    cores <- 1
  }
  
  ret <- fit.robmixglm(X,Y,family,offset=offset,gh=norm.gauss.hermite(quadpoints),notrials,EMTol,cores,verbose)

  class(ret) <- "robmixglm"

  ret$call <- call
  ret$family <- family
  
  ret$X <- X
  ret$Y <- Y
  ret$offset <- offset
  
  ret$model <- mf
  ret$terms <- mt

  ret$xlevels <- .getXlevels(mt, mf)

  ret$quadpoints <- quadpoints
  ret$notrials <- notrials
  ret$EMTol <- EMTol
  ret$verbose <- verbose
  return(ret)
}

