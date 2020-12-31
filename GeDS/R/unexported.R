getbases <- function(obj, n = 2L){
  n <- as.integer(n)
  if(!(n %in% 2L:4L) ) stop("Please set the order between 2 and 4") else str <- NULL
  if(n==2L) str <- "Linear"
  if(n==3L) str <- "Quadratic"
  if(n==4L) str <- "Cubic"
  bases <- eval(parse(text=paste0("obj$",str,"$Basis")))
  return(bases)
}

buildbasis <- function(parm, object, n){
  kn <- knots(object,n=n, options = "all")
  mat <- splineDesign(knots = kn, x=parm, ord=n)
  return(mat)
}

confint.GeDS <- function(object, parm, level = 0.95, n = 3L, ...){
  n <- as.integer(n)
  if(!(n %in% 2L:4L) ) stop("Please set the order between 2 and 4") else str <- NULL
  flag <- TRUE
  env <- object$Args
  extr <- env$extr
  Xmat <- getbases(object, n=n)
  m <- NROW(Xmat)
  K <- NCOL(Xmat)
  if(n==2L) str <- "Linear"
  if(n==3L) str <- "Quadratic"
  if(n==4L) str <- "Cubic"
  pred <- if(missing(parm)) eval(parse(text=paste0("object$",str,"$Predicted"))) else
    predict.GeDS(object,newdata=data.frame(parm),type="link")
  X <- if(missing(parm)) Xmat else buildbasis(parm,object,n)

  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3) # stats package


  if(object$Type=="GLM - Univ"){
    flag <- FALSE
    family <- env$family
    norm <- function(x,extr){
      x <- (x-extr[1])/diff(extr)
    }
    #variance
    alpha_hat <- coef(object, n=n, onlySpline=T)
    if(!is.canonical.link(family))
      stop("Confidence intervals can be computed only if the canonical link is used")
    if(family$family %in% c("poisson","quasipoisson")){
      c2 <- exp
    }
    if(family$family =="Gamma"){
      c2 <- function(x) -(1/x^2)
    }
    if(family$family == "binomial"){
      c2 <- function(x) exp(x)/(1+exp(x))^2
    }
    if(family$family %in% c("gaussian","quasi")){
      c2 <- function(x) 1
    }
    if(family$family == "inverse.gaussian"){
      c2 <- function(x) -2/x^3
    }

    What <- diag(as.numeric(c2(Xmat%*%alpha_hat)))
    Gmat <- (m^-1)*t(Xmat)%*%What%*%Xmat
    Ginv <- solve(Gmat)
    variance <- (K-n)^-1*X%*%Ginv%*%t(X)

    if(length(variance>1)) variance <- diag(variance)
    #bias
    # we need one more degree in order to estimate the derivative
    if(n==4L){
      kn <- makenewknots(knots(object,n=2L, options="internal"),5)
      env$kn <- kn
      env$predicted <- object$Cubic$Predicted
      must <- family$linkinv(env$predicted);
      reg <- SplineReg_GLM(X=env$X,Y=env$Y,Z=env$Z,offset=env$offset,
                           weights=env$weights,extr=env$extr,InterKnots=env$kn,n=5L,
                           family=env$family,mustart=must)
      # reg <- evalq({
      #   must <- family$linkinv(predicted);
      #   SplineReg_GLM(X=X,Y=Y,Z=Z,offset=offset,
      #                 weights=weights,extr=extr,InterKnots=kn,n=5L,
      #                 family=family,mustart=must)
      # },
      # envir = env)
      thetas <- reg$Theta
    } else {
      thetas <- coef(object, n=n+1L)
      kn <- knots(object,n=n+1L,options="internal")
    }
    #  derive
    #Bernoulli
    kn <- norm(kn,extr)

    if(missing(parm)) parm <- env$X
    parm2 <- norm(parm,extr)

    basis <- splines::splineDesign(knots=sort(c(kn, rep(env$extr,n+1))), x=parm2, ord = n+1, derivs = n, outer.ok = TRUE)
    #probably better to use even higher degree,
    der <- as.numeric(basis%*%thetas)
    kn2 <- sort(c(kn,0,1))

    bias <- numeric(length(parm2))
    par<- NULL

    for(i in 1:length(parm2)){
      par <- parm2[i]
      id <- max(which(kn2<par))
      bernoulli <- BerPoly(n+1,(par-kn2[id])*length(kn)) # sum
      bias[i]<- -der[i]/(factorial(n+1)*length(kn)^(n+1))*bernoulli
    }

    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    #print(bias)
    ci[] <- as.numeric(pred - bias) + sqrt(variance) %o% qnorm(a)
  }
  if(object$Type=="LM - Univ"){
    flag <- FALSE
    matcb <- matrix(0,K,K)
    for(i in 1:m){
      matcb <- matcb + Xmat[i,]%*%t(Xmat[i,])
    }
    matcb <- matcb/m
    matcbinv <- solve(matcb)

    reg <- eval(parse(text=paste0("object$",str)))
    df <- reg$temporary$df
    resid <- reg$Residuals
    sigma_hat <- sqrt(sum(resid^2)/df)

    if(missing(parm)) parm <- env$X
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    ci[] <- as.numeric(pred) + m^(-.5)*diag(sigma_hat*Xmat%*%matcbinv%*%t(Xmat)) %o% qnorm(a)
  }

  if(flag) stop("Function available only in the unviariate framework")

  return(ci)
}

is.canonical.link <- function(family){
  str <- family$family
  fam <- eval(parse(text=paste0(str,"()")))
  res <- fam$link==family$link
  return(res)
}

BerPoly <- function(n,x){
  v <- 0:n
  res <- choose(n,v) * Bernoulli(n-v) * x^v
  res <- sum(as.numeric(res))
  return(res)
}

# just copy and paste from package::stats
# unexported there, same here
format.perc <- function(probs, digits)
  ## Not yet exported, maybe useful in other contexts:
  ## quantile.default() sometimes uses a version of it
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
        "%")

gedsint <- function(val, knts, coefs, n){
  if(val<=min(knts)) return(0)
  if(val>max(knts)) return(gedsint(max(knts), knts=knts, coefs=coefs, n))
  pos <- min(which(knts>=val))
  matrice <- splines::splineDesign(knots=c(knts,max(knts)),
                                   derivs=rep(0,length(val)),x=val, ord=n+1, outer.ok = T)

  ris <- as.numeric(matrice[,1:(pos-1)]%*%coefs[1:(pos-1)])
  return(ris)
}


