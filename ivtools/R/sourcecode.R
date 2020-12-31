globalVariables(c(".", ".SD", "p00.0", "p00.1", "p00.2",
  "p01.0", "p01.1", "p01.2",
  "p10.0", "p10.1", "p10.2",
  "p11.0", "p11.1", "p11.2"))
  
aggr <- function(x, clusters){
  #seems like lapply(.SD, sum) on data.table drops 
  #columns with identical names, thus ensure that no identical names
  colnames(x) <- NULL 
  temp <- data.table(x)
  temp <- as.matrix(temp[, j=lapply(.SD, sum), by=clusters])[, -1]

} 

ah <- function(formula, data, weights, robust=FALSE){ 

  if(missing(weights))
    weights <- rep(1, nrow(data))
  mf <- model.frame(formula=formula, data=data)
  incl.names <- rownames(mf)
  incl.rows <- match(incl.names, rownames(data))  
  data <- data[incl.rows, ]
  weights <- weights[incl.rows]
  surv <- mf[, 1]
  X <- model.matrix(object=formula, data=data)[, -1, drop=FALSE]
  fit <- ahaz(surv=surv, X=X, weights=weights, robust=robust)
  fit$call <- match.call()
  fit$formula <- eval(formula)
  coef <- solve(fit$D, fit$d)
  coeff.names <- colnames(X)
  names(coef) <- coeff.names
  fit$coefficients <- coef
  vcov <- solve(fit$D)%*%fit$B%*%solve(fit$D)
  rownames(vcov) <- coeff.names
  colnames(vcov) <- coeff.names
  fit$vcov <- vcov
  fit$incl <- incl.names
  class(fit) <- c("ah", "ahaz")
  return(fit) 
   
}

confint.ivmod <- function(object, parm, level=0.95, ...){
  
  est <- object$est
  se <- sqrt(diag(object$vcov))
  qqq <- abs(qnorm((1-level)/2))
  lower <- est-qqq*se
  upper <- est+qqq*se
  ci <- cbind(lower, upper)
  colnames(ci) <- paste(100*c((1-level)/2 , (1+level)/2), "%")
  return(ci)
  
}

estfun <- function(object, lower, upper, step){
 
  #extract stuff from fitted gest object 
  input <- object$input
  data <- input$data
  n <- nrow(data)
  fitZ.L <- input$fitZ.L
  fitX.LZ <- input$fitX.LZ
  fitX.L <- input$fitX.L
  if(is.null(fitX.LZ)){
    Z <- all.vars(fitZ.L$formula)[1]
    X <- input$X     
  }else{
    Z <- input$Z
    X <- all.vars(fitX.LZ$formula)[1]    
  }
  Y <- input$Y
  if(inherits(x=object, what="ivglm")){
    type <- input$link
    fitY.LZX <- input$fitY.LZX
  }
  if(inherits(x=object, what="ivcoxph")){
    type <- "coxph"
    fitY.LZX <- input$fitT.LZX
    fit.detail <- coxph.detail(object=fitY.LZX)
  }
  if(inherits(x=object, what="ivah"))
    stop("estfun is not implemented for ivah objects.")
  formula <- input$formula
  y <- object$t
   
  #collect useful stuff  
  stuff <- utilityfun(Z=Z, X=X, Y=Y, type=type, data=data, 
    formula=formula, y=y, fitY.LZX=fitY.LZX, fitX.LZ=fitX.LZ, fitX.L=fitX.L, 
    fitZ.L=fitZ.L, fit.detail=fit.detail)
  est.nuisance <- stuff$est.nuisance
  nZ <- stuff$nZ
  nX.LZ <- stuff$nX.LZ
  nX.L <- stuff$nX.L 
  nY <- stuff$nY
  npsi <- stuff$npsi
  designZ <- stuff$designZ
  designX.LZ <- stuff$designX.LZ
  designX.L <- stuff$designX.L
  designY <- stuff$designY
  designpsi <- stuff$designpsi
  weights <- stuff$weights
  
  if(object$converged)  
    est <- object$est
  else
    est <- rep(0, length(object$est))
  if(missing(lower))
    lower <- est-0.5
  if(missing(upper))
    upper <- est+0.5
  if(missing(step))
    step <- rep(0.01, npsi)
    
  f <- list()
  for(i in 1:npsi){
    tmp <- seq(lower[i], upper[i], step[i])
    m <- matrix(nrow=length(tmp), ncol=2)
    colnames(m) <- c("psi", "Hsum")
    m[, 1] <- tmp 
    for(j in 1:length(tmp)){
      psii <- est
      psii[i] <- tmp[j]
      Uj <- Upsifun(b=c(psii, est.nuisance), Z=Z, X=X, Y=Y, type=type, 
        fitZ.L=fitZ.L, fitX.LZ=fitX.LZ, fitX.L=fitX.L, fitY.LZX=fitY.LZX, 
        npsi=npsi, nZ=nZ, nX.LZ=nX.LZ, nX.L=nX.L, nY=nY, designpsi=designpsi, 
        designZ=designZ, designX.LZ=designX.LZ, designX.L=designX.L, 
        designY=designY, weights=weights, data=data)
      m[j, 2] <- colSums(Uj)[i] 
    }
    f[[i]] <- m
  }
  names(f) <- names(est)
  
  out <- list(f=f, est=est)
  class(out) <- "estfun"
  return(out) 
} 

expand <- function(x, names){
  if(is.vector(x))
    x <- x[match(names, names(x))]
  if(is.matrix(x))
    x <- x[match(names, rownames(x)), , drop=FALSE]
  return(x)  
} 

gest <- function(Z, X, Y, type, data, formula=~1, y=NULL, fitY.LZX=NULL, 
  fitX.LZ=NULL, fitX.L=NULL, fitZ.L=NULL, clusterid=NULL, vcov.fit=vcov.fit, 
  fit.detail=NULL, ...){
  
  dots <- list(...)
  n <- nrow(data) 
  
  #find optimal y if not provided
  if(type=="coxph" & is.null(y)){
    f <- function(y){
      fit <- gest(Z=Z, X=X, Y=Y, type="coxph", formula=formula, y=y, 
        fitY.LZX=fitY.LZX, fitX.LZ=fitX.LZ, fitX.L=fitX.L, fitZ.L=fitZ.L, 
        data=data, clusterid=clusterid, vcov.fit=vcov.fit, 
        fit.detail=fit.detail, ...)
      return(sum(diag(fit$vcov)))
    }
    if(inherits(x=fitY.LZX, what="survfit")){
      ss <- summary(fitY.LZX)
      time <- ss$time
    }
    if(inherits(x=fitY.LZX, what="coxph")){
      time <- fit.detail$time
    }
    interval <- c(min(time), max(time))
    y <- optimize(f=f, interval=interval)$minimum  
  }
 
  #number of clusters?
  if(!is.null(clusterid))
    ncluster <- length(unique(data[, clusterid]))   
  
  #collect useful stuff  
  stuff <- utilityfun(Z=Z, X=X, Y=Y, type=type, data=data, 
    formula=formula, y=y, fitY.LZX=fitY.LZX, fitX.LZ=fitX.LZ, fitX.L=fitX.L,
    fitZ.L=fitZ.L, fit.detail=fit.detail)
  est.nuisance <- stuff$est.nuisance
  nZ <- stuff$nZ
  nX.LZ <- stuff$nX.LZ
  nX.L <- stuff$nX.L
  nY <- stuff$nY
  npsi <- stuff$npsi
  designZ <- stuff$designZ
  designX.LZ <- stuff$designX.LZ
  designX.L <- stuff$designX.L
  designY <- stuff$designY
  designpsi <- stuff$designpsi
  weights <- stuff$weights 
  t1 <- stuff$t1 
  
  #help function for summing estimating functions
  estfunsums <- function(b){
    if(length(b)==npsi)
      b <- c(b, est.nuisance)
    return(colSums(Upsifun(b=b, Z=Z, X=X, Y=Y, type=type, 
      fitZ.L=fitZ.L, fitX.LZ=fitX.LZ, fitX.L=fitX.L, fitY.LZX=fitY.LZX, 
      npsi=npsi, nZ=nZ, nX.LZ=nX.LZ, nX.L=nX.L, nY=nY, designpsi=designpsi, 
      designZ=designZ, designX.LZ=designX.LZ, designX.L=designX.L, 
      designY=designY, weights=weights, data=data)))
  }
 
  #compute estimate of psi
  args <- list(fn=estfunsums)
  args[names(dots)] <- dots
  if(is.na(match("x", names(args))))
    args$x <- rep(0, npsi)  
  args$global <- "cline" 
  fit.nleqslv <- do.call("nleqslv", args=args)
  est <- fit.nleqslv$x
  converged <- fit.nleqslv$termcd==1
  
  #if psi is scalar and no solution to estimating equation was found,
  #then try a range of different starting values around 0
  if(npsi==1 & !converged){
    if(args$x==0) 
      x <- -1
    else
      x <- 0
    xmax <- 10
    while(!converged & abs(x)<=xmax){
      args$x <- x
      fit.nleqslv <- do.call("nleqslv", args=args)
      est <- fit.nleqslv$x
      converged <- fit.nleqslv$termcd==1
      if(x<0)
        x <- -x
      else
        x <- -x-1  
    }
  }
  if(!converged)
    est <- rep(NA, npsi) 
  names(est) <- paste0(X, ":", colnames(designpsi))
  names(est)[1] <- X
   
 #compute variance if requested, and if solution to estimating equation was found
  if(vcov.fit & converged){   
    if(!is.null(fitZ.L))
      sandwich.fitZ.L <- sandwich(fit=fitZ.L, data=data, weights=weights)
    if(!is.null(fitX.LZ))
      sandwich.fitX.LZ <- sandwich(fit=fitX.LZ, data=data, weights=weights) 
    if(!is.null(fitX.L))
      sandwich.fitX.L <- sandwich(fit=fitX.L, data=data, weights=weights) 
    if(!is.null(fitY.LZX))
      sandwich.fitY.LZX <- sandwich(fit=fitY.LZX, data=data, weights=weights,
        t=y, fit.detail=fit.detail) 
    Upsi <- Upsifun(b=c(est, est.nuisance), Z=Z, X=X, Y=Y, type=type, 
      fitZ.L=fitZ.L, fitX.LZ=fitX.LZ, fitX.L=fitX.L, fitY.LZX=fitY.LZX, 
      npsi=npsi, nZ=nZ, nX.LZ=nX.LZ, nX.L=nX.L, nY=nY, designpsi=designpsi, 
      designZ=designZ, designX.LZ=designX.LZ, designX.L=designX.L, 
      designY=designY, weights=weights, data=data)
    Ipsi <- jacobian(func=estfunsums, x=c(est, est.nuisance))/n
    if(is.null(fitX.LZ)){
      Ud <- sandwich.fitZ.L$U
      Id <- cbind(matrix(0, nrow=nZ, ncol=npsi),
        sandwich.fitZ.L$I,
        matrix(0, nrow=nZ, ncol=nY))
      nd <- nZ
      namesd <- names(fitZ.L$coefficients)  
    }else{
      Ud <- cbind(sandwich.fitX.LZ$U, sandwich.fitX.L$U) 
      Id <- rbind(cbind(matrix(0, nrow=nX.LZ, ncol=npsi),
        sandwich.fitX.LZ$I,
        matrix(0, nrow=nX.LZ, ncol=nX.L),
        matrix(0, nrow=nX.LZ, ncol=nY)),
        cbind(matrix(0, nrow=nX.L, ncol=npsi),
        matrix(0, nrow=nX.L, ncol=nX.LZ),
        sandwich.fitX.L$I,
        matrix(0, nrow=nX.L, ncol=nY)))
      nd <- nX.LZ+nX.L
      namesd <- c(names(fitX.LZ$coefficients), names(fitX.L$coefficients))  
    }  
    if(type=="identity" | type=="log"){
      UY <- NULL
      IY <- NULL
    }
    if(type=="logit" | type=="coxph"){
      UY <- sandwich.fitY.LZX$U
      IY <- cbind(matrix(0, nrow=nY, ncol=npsi+nd), sandwich.fitY.LZX$I)
    }
    U <- cbind(Upsi, Ud, UY)
    I <- rbind(Ipsi, Id, IY)
    if(is.null(clusterid)){
      J <- var(U, na.rm=TRUE)
      vcov <- (solve(I)%*%J%*%t(solve(I))/n)
    }
    else{
      U <- aggr(x=U, clusters=data[, clusterid])
      J <- var(U, na.rm=TRUE)
      vcov <- (solve(I)%*%J%*%t(solve(I))*ncluster/n^2)
    }   
    vcov <- as.matrix(vcov[1:npsi, 1:npsi])
  }
  else{
    U <- NA
    I <- NA
    vcov <- matrix(NA, nrow=npsi, ncol=npsi)  
  }
  rownames(vcov) <- names(est)
  colnames(vcov) <- names(est)
  
  result <- list(est=est, vcov=vcov, estfunall=U, d.estfun=I, 
    converged=converged, fitZ.L=fitZ.L, fitX.L=fitX.L, y=y)
  
  class(result) <- "gest"

  return(result)

}

Hfun <- function(fit, data, fit.detail){
  if(inherits(x=fit, what="survfit")){
    #need to use summary(fit), since n.events and n.risk from fit 
    #behave strange when there is left truncation     
    ss <- summary(fit)
    strata <- ss$strata
    #n.strata is unweighted
    n.strata <- summary(strata)
    K <- length(n.strata)
    names.strata <- names(n.strata)
    time <- ss$time
    #n.event and n.risk are weighted
    n.event <- ss$n.event
    n.risk <- ss$n.risk  
    dH <- n.event/n.risk
    H <- list()
    breaks <- c(0, cumsum(n.strata))
    for(k in 1:K){
      incl <- (breaks[k]+1):breaks[k+1] 
      H[[k]] <- stepfun(time[incl], c(0, cumsum(dH[incl])))
    }
    names(H) <- names.strata
  }
  if(inherits(x=fit, what="coxph")){  
    time <- fit.detail$time
    #dH is weighted
    dH <- fit.detail$hazard
    H <- stepfun(time, c(0, cumsum(dH)))         
  }   
  
  return(H)
}


ivah <- function(estmethod, X, T, fitZ.L=NULL, fitX.LZ=NULL, fitT.LX=NULL,
  data, ctrl=FALSE, clusterid=NULL, event, max.time, max.time.psi, n.sim=100, 
  vcov.fit=TRUE, ...){
  
  call <- match.call()
  input <- as.list(environment())
  
  if(estmethod=="ts"){
    if(is.null(fitT.LX))
      stop("TS estimation requires fitT.LX") 
    if(is.null(fitX.LZ))
      stop("TS estimation requires fitX.LZ")
    if(!is.null(fitZ.L))
      warning("TS estimation doesn't use fitZ.L, this will be ignored")
    result <- tsest(fitX.LZ=fitX.LZ, fitY.LX=fitT.LX, data=data, ctrl=ctrl, 
      clusterid=clusterid, vcov.fit=vcov.fit)
  }
  if(estmethod=="g"){
    if(is.null(fitZ.L))
      stop("G-estimation requires fitZ.L")
    if(missing(X))
      stop("G-estimation requires the name of the exposure X")
    if(missing(T))
      stop("G-estimation requires the name of the follow-up time T")
    if(missing(max.time))
      max.time <- max(data[, T])
    if(missing(max.time.psi))
      max.time.psi <- max(data[, T])
    result <- scs(X=X, time=T, status=event, data=data, fitZ.L=fitZ.L, 
      max.time=max.time, max.time.beta=max.time.psi, n.sim=n.sim, ...)
    result$converged <- !is.na(result$est)
    
  }
  if(!result$converged)
    warning("No solution to the estimating equation was found")
  result$call <- call 
  result$input <- input  
  class(result) <- c("ivah", "ivmod")
  return(result)

}

ivbounds <- function(data, Z, X, Y, monotonicity=FALSE, weights){
  
  call <- match.call()
  
  zlev <- c(0, 1)
  xlev <- c(0, 1)
  ylev <- c(0, 1) 
  if(is.data.frame(data)){
    if(missing(weights))
      weights <- rep(1, nrow(data))
    incl <- !is.na(data[, Z]) & !is.na(data[, X]) & !is.na(data[, Y])
    data <- data[incl, ]
    weights <- weights[incl]
    if(length(unique(data[, Z]))==3)
      zlev <- c(zlev , 2)
  }
  else{
    if(missing(weights))
      weights <- rep(1, )
    if(length(data)==12)
      zlev <- c(zlev , 2)
  }
  for(y in ylev){
    for(x in xlev){
      for(z in zlev){
        p <- paste0("p", y, x, ".", z)
        if(is.data.frame(data)){
          num <- weighted.mean(x=(data[, Y]==y & data[, X]==x & data[, Z]==z),
            w=weights)
          den <- weighted.mean(x=(data[, Z]==z), w=weights)
          val <- num/den
        }
        else{
          val <- data[p]
        }
        assign(x=p, value=val)
      }
    }
  }   
  if(length(zlev)==2){
    if(monotonicity){
      p0.low <- "p10.0"
      p0.upp <- "1-p00.0"
      p1.low <- "p11.1"
      p1.upp <- "1-p01.1"
    }else{
      p0.low <- c("p10.0+p11.0-p00.1-p11.1", 
        "p10.1", 
        "p10.0", 
        "p01.0+p10.0-p00.1-p01.1")
      p0.upp <- c("p01.0+p10.0+p10.1+p11.1", 
        "1-p00.1", 
        "1-p00.0", 
        "p10.0+p11.0+p01.1+p10.1") 
      p1.low <- c("p11.0", 
        "p11.1", 
        "-p00.0-p01.0+p00.1+p11.1", 
        "-p01.0-p10.0+p10.1+p11.1")
      p1.upp <- c("1-p01.1", 
        "1-p01.0", 
        "p00.0+p11.0+p10.1+p11.1", 
        "p10.0+p11.0+p00.1+p11.1")   
    }
      
  }
  if(length(zlev)==3){
    if(monotonicity){
      p0.low <- "p10.0"
      p0.upp <- "1-p00.0"
      p1.low <- "p11.2"
      p1.upp <- "1-p01.2"  
    }else{
      p0.low <- c("p10.0", 
        "p10.1", 
        "p10.2", 
        "p10.0+p11.0+p10.1+p01.1-1", 
        "p10.0+p01.0+p10.1+p11.1-1", 
        "p10.1+p11.1+p10.2+p01.2-1", 
        "p10.1+p01.1+p10.2+p11.2-1", 
        "p10.2+p11.2+p10.0+p01.0-1", 
        "p10.2+p01.2+p10.0+p11.0-1")
      p0.upp <- c("1-p00.0", 
        "1-p00.1", 
        "1-p00.2", 
        "p10.0+p01.0+p10.1+p11.1", 
        "p10.0+p11.0+p10.1+p01.1", 
        "p10.1+p01.1+p10.2+p11.2", 
        "p10.1+p11.1+p10.2+p01.2", 
        "p10.2+p01.2+p10.0+p11.0", 
        "p10.2+p11.2+p10.0+p01.0")
      p1.low <- c("p11.0", 
        "p11.1", 
        "p11.2", 
        "p10.0+p11.0-p10.1-p01.1", 
        "-p10.0-p01.0+p10.1+p11.1", 
        "p10.1+p11.1-p10.2-p01.2", 
        "-p10.1-p01.1+p10.2+p11.2", 
        "p10.2+p11.2-p10.0 -p01.0", 
        "-p10.2-p01.2+p10.0+p11.0")
      p1.upp <- c("1-p01.0", 
        "1-p01.1", 
        "1-p01.2", 
        "p10.0+p11.0-p10.1-p01.1+1", 
        "-p10.0-p01.0+p10.1+p11.1+1", 
        "p10.1+p11.1-p10.2-p01.2+1", 
        "-p10.1-p01.1+p10.2+p11.2+1", 
        "p10.2+p11.2-p10.0-p01.0+1", 
        "-p10.2-p01.2+p10.0+p11.0+1")
    }   
  }

  p0.low.num <- sapply(p0.low, function(x) eval(parse(text=x))) 
  p0.upp.num <- sapply(p0.upp, function(x) eval(parse(text=x)))
  p1.low.num <- sapply(p1.low, function(x) eval(parse(text=x)))
  p1.upp.num <- sapply(p1.upp, function(x) eval(parse(text=x)))
  p0 <- c(max(p0.low.num), min(p0.upp.num))
  p1 <- c(max(p1.low.num), min(p1.upp.num))  
  names(p0) <- (c("min", "max"))
  names(p1) <- (c("min", "max"))
  p0.symbolic <- c(p0.low[which.max(p0.low.num)], p0.upp[which.min(p0.upp.num)])
  p1.symbolic <- c(p1.low[which.max(p1.low.num)], p1.upp[which.min(p1.upp.num)])
  names(p0.symbolic) <- (c("min", "max"))
  names(p1.symbolic) <- (c("min", "max"))
  
  IVinequality <- TRUE
  conditions <- NULL
  if(length(zlev)==2)
    temp <- cbind(c(0, 0, 1, 1), c(0, 1, 0, 1))
  if(length(zlev)==3)
    temp <- cbind(c(0, 0, 0, 1, 1, 1), c(0, 1, 2, 0, 1, 2))
  for(x in xlev){
    for(i in 1:nrow(temp)){
      cond <- paste0("p", 0, x, ".", temp[i, 1], "+",
        "p", 1, x, ".", temp[i, 2])  
      cond <- paste0(cond, "<=1")
      if(!eval(parse(text=cond))){
        IVinequality <- FALSE
        conditions <- c(conditions, cond)
      }
    }     
  }
  result <- list(call=call, p0=p0, p1=p1, p0.symbolic=p0.symbolic, 
    p1.symbolic=p1.symbolic, IVinequality=IVinequality, conditions=conditions)
  class(result) <- "ivbounds"
    
  return(result)

}

ivcoxph <- function(estmethod, X, fitZ.L=NULL, fitX.LZ=NULL, fitX.L=NULL, 
  fitT.LX=NULL, fitT.LZX=NULL, data, formula=~1, ctrl=FALSE, clusterid=NULL, 
  t=NULL, vcov.fit=TRUE, ...){
  
  call <- match.call()
  input <- as.list(environment())
  
  if(estmethod=="ts"){
    if(is.null(fitT.LX))
      stop("TS estimation requires fitT.LX") 
    if(is.null(fitX.LZ))
      stop("TS estimation requires fitX.LZ")
    if(!is.null(fitZ.L))
      warning("TS estimation doesn't use fitZ.L, this will be ignored")
    if(!is.null(fitX.L))
      warning("TS estimation doesn't use fitX.L, this will be ignored") 
    if(!is.null(fitT.LZX))
      warning("TS estimation doesn't use fitT.LZX, this will be ignored")
    result <- tsest(fitX.LZ=fitX.LZ, fitY.LX=fitT.LX, data=data, ctrl=ctrl, 
      clusterid=clusterid, vcov.fit=vcov.fit)
  }
  if(estmethod=="g"){
    if(is.null(fitZ.L) & is.null(fitX.LZ))
      stop("G-estimation requires either fitZ.L or fitX.LZ")
    if(!is.null(fitX.LZ) & is.null(fitX.L))
      stop("G-estimation with fitX.LZ requires fitX.L as well")
    if(!is.null(fitZ.L) & !is.null(fitX.LZ) & !is.null(fitX.L))
      warning("G-estimation doesn't require both fitZ.L, fitX.LZ and fitX.L, 
        only fitX.LZ and fitX.L will be used")
    if(!is.null(fitT.LX))
      warning("G-estimation doesn't use fitT.LX, this will be ignored")
    if(is.null(fitT.LZX))
      stop("G-estimation requires fitT.LZX")
    if(is.null(fitX.LZ)){
      Z <- all.vars(fitZ.L$formula)[1] 
      if(missing(X))
        stop("G-estimation with fitZ.L requires the name of the exposure X")   
    }else{
      Z <- NULL
      X <- all.vars(fitX.LZ$formula)[1]
    } 
    type <- "coxph"
    if(inherits(x=fitT.LZX, what="coxph")){
      tmp <- all.vars(fitT.LZX$formula[[2]])
      fit.detail <- coxph.detail(object=fitT.LZX)
    }
    if(inherits(x=fitT.LZX, what="survfit")){
      #survfit object has no formula element, so need to get it from call,
      #need to use eval, since the fit$call$formula will be literary what the user
      #gave as argument, e.g. if formula=f, then fit$call$formula is f, not the 
      #formula contained in f
      tmp <- all.vars(eval(fitT.LZX$call$formula)[[2]]) 
      fit.detail <- NULL  
    }
    Y <- tmp[length(tmp)-1]
    result <- gest(Z=Z, X=X, Y=Y, type=type, data=data, formula=formula, y=t,
      fitY.LZX=fitT.LZX, fitX.LZ=fitX.LZ, fitX.L=fitX.L, fitZ.L=fitZ.L, 
      clusterid=clusterid, vcov.fit=vcov.fit, fit.detail=fit.detail, ...)
    result$t <- result$y
    result$y <- NULL
  }
  if(!result$converged)
    warning("No solution to the estimating equation was found")
  result$call <- call 
  result$input <- input
  class(result) <- c("ivcoxph", "ivmod")
  return(result)

}

ivglm <- function(estmethod, X, Y, fitZ.L=NULL, fitX.LZ=NULL, fitX.L=NULL, 
  fitY.LX=NULL, fitY.LZX=NULL, data, formula=~1, ctrl=FALSE, clusterid=NULL, 
  link, vcov.fit=TRUE, ...){
   
  call <- match.call()
  input <- as.list(environment())
  
  if(estmethod=="ts"){
    if(is.null(fitY.LX))
      stop("TS estimation requires fitY.LX") 
    if(is.null(fitX.LZ))
      stop("TS estimation requires fitX.LZ")
    if(!is.null(fitZ.L))
      warning("TS estimation doesn't use fitZ.L, this will be ignored")
    if(!is.null(fitX.L))
      warning("TS estimation doesn't use fitX.L, this will be ignored") 
    if(!is.null(fitY.LZX))
      warning("TS estimation doesn't use fitY.LZX, this will be ignored")
    result <- tsest(fitX.LZ=fitX.LZ, fitY.LX=fitY.LX, data=data, ctrl=ctrl, 
      clusterid=clusterid, vcov.fit=vcov.fit)
  }
  if(estmethod=="g"){
    if(is.null(fitZ.L) & is.null(fitX.LZ))
      stop("G-estimation requires either fitZ.L or fitX.LZ")
    if(!is.null(fitX.LZ) & is.null(fitX.L))
      stop("G-estimation with fitX.LZ requires fitX.L as well")
    if(!is.null(fitZ.L) & !is.null(fitX.LZ) & !is.null(fitX.L))
      warning("G-estimation doesn't require both fitZ.L, fitX.LZ and fitX.L, 
        only fitX.LZ and fitX.L will be used")
    if(!is.null(fitY.LX))
      warning("G-estimation doesn't use fitY.LX, this will be ignored")
    if(is.null(fitY.LZX) & link=="logit")
      stop("G-estimation with logit link requires fitY.LZX") 
    if(is.null(fitX.LZ)){
      Z <- all.vars(fitZ.L$formula)[1] 
      if(missing(X))
        stop("G-estimation with fitZ.L requires the name of the exposure X")   
    }else{
      Z <- NULL
      X <- all.vars(fitX.LZ$formula)[1]
    }    
    if(missing(Y) & link=="identity")
      stop("G-estimation with identity link requires the name of the outcome Y")
    if(missing(Y) & link=="log")
      stop("G-estimation with log link requires the name of the outcome Y")
    if(link=="logit")
      Y <- all.vars(fitY.LZX$formula)[1]
    type <- link
    result <- gest(Z=Z, X=X, Y=Y, type=type, data=data, formula=formula,
      fitY.LZX=fitY.LZX, fitX.LZ=fitX.LZ, fitX.L=fitX.L, fitZ.L=fitZ.L, 
      clusterid=clusterid, vcov.fit=vcov.fit, ...)
  }
  if(!result$converged)
    warning("No solution to the estimating equation was found")
  result$call <- call 
  result$input <- input
  class(result) <- c("ivglm", "ivmod")
  return(result)

}

plot.estfun <- function(x, ...){
  dots <- list(...)
  args <- list(type="l")
  args[names(dots)] <- dots
  if(is.na(match("xlab", names(args))))
    args$xlab <- expression(psi)
  if(is.na(match("ylab", names(args))))
    args$ylab <- expression(H(psi)) 
  npar <- length(x$est)
  nr <- ceiling(sqrt(npar))
  nc <- ceiling(npar/nr)
  par(mfrow=c(nr, nc), mar=c(4.2,5.5,2,0.5))
  for(i in 1:npar){
    args$x <- x$f[[i]][, 1]
    args$y <- x$f[[i]][, 2] 
    args$main <- names(x$est)[i]
    do.call("plot", args=args)
    abline(h=0, lty="dashed")
    abline(v=x$est[i], lty="dashed")
  }
}

plot.ivah <- function (x, gof=FALSE, CI.level=0.95, ...){
  
  dots <- list(...)
  if(gof==TRUE){
    res <- x$GOF.resamp
    ylim <- c(min(c(res[2, ], res[2:22, ])), max(c(res[2, ], res[2:22, ])))
    args <- list(x=res[1, ], y=res[2, ], type="n", ylim=ylim)
    args[names(dots)] <- dots
    if(is.na(match("xlab", names(args))))
      args["xlab"] <- "Time"
    if(is.na(match("ylab", names(args))))
      args["ylab"] <- "Goodness of fit process"
    do.call("plot", args=args)   
    matlines(t(res[1, , drop=FALSE]), t(res[2:20, ]), col="grey", lty="solid")
    lines(res[1, ], res[2, ])
    #abline(0, 0)	
  }
  else{
    qqq <- abs(qnorm((1-CI.level)/2))
    B <- x$B
    se <- x$se_B 
    l <- B-qqq*se
    u <- B+qqq*se
    beta <- x$est
    stime <- x$stime
    ylim <- c(min(c(l, beta*stime)), max(c(u, beta*stime)))
    args <- list(x=stime, y=B, type="l", ylim=ylim)
    args[names(dots)] <- dots
    if(is.na(match("xlab", names(args))))
      args["xlab"] <- "Time"
    if(is.na(match("ylab", names(args))))
      args$ylab <- expression(B(t)) 
    do.call("plot", args=args)
    lines(stime, u, lty="dashed")
    lines(stime, l, lty="dashed")
    #abline(0, 0)
    lines(stime, beta*stime)
  }
}


print.ivmod <- function (x, digits=max(3L, getOption("digits")-3L), ...) 
{
    if (length(x$est)) {
        cat("\nCoefficients:\n")
        print.default(format(x$est, digits=digits), print.gap=2, 
            quote=FALSE)
        cat("\n")
    }
    else {
        cat("No coefficients\n\n")
    }
}

print.summary.ivbounds <- function(x, digits=max(3L, getOption("digits")-3L), 
   ...) {
  cat("\nCall:  ", "\n", paste(deparse(x$call), sep="\n", collapse="\n"), 
      "\n", sep="")
  if(x$IVinequality){
    cat("\nThe IV inequality is not violated\n\n")
  }else{
    cat("\nThe IV inequality is violated at the following conditions:\n")
    for(i in length(x$conditions)){
      print(x$conditions[i], quote=FALSE)
      cat("  \n")
    }
  }
  cat("Symbolic bounds:\n")
  print(x$coefficients.symbolic, quote=FALSE)
  cat("\n")
  cat("Numeric bounds:\n")
  print(x$coefficients, digits)
  cat("\n") 
}

print.summary.ivmod <- function(x, digits=max(3L, getOption("digits")-3L), 
  signif.stars=getOption("show.signif.stars"), ...) {
  cat("\nCall:  ", "\n", paste(deparse(x$call), sep="\n", collapse="\n"), 
      "\n", sep="")
  if(!is.null(x$t))
    cat("\nEquation solved at t =", x$t, "\n")
  if(!is.null(x$test0)){
    cat("\nTest for non-significant exposure effect. H_0: B(t)=0 \n")
    cat("   \n")
    prmatrix(signif(x$test0, digits))
    cat("   \n")
    cat("Goodness-of-fit test for constant effects model\n")
    prmatrix(signif(x$test_gof, digits))
    cat("   \n");  
    cat("Constant effect model  \n")
  } 
  cat("\nCoefficients:", "\n")
  printCoefmat(x$coefficients, digits=digits, signif.stars=signif.stars, 
    na.print = "NA", ...)
  cat("\n") 
}

sandwich <- function(fit, data, weights, t, fit.detail){

  n <- nrow(data)
  if(missing(weights))
    weights <- rep(1, n)
  
  if(inherits(x=fit, what="glm")){
    
    #---meat---
    
    m <- expand(model.matrix(fit), rownames(data))
    res <- expand(residuals(fit, type="response"), rownames(data)) 
    U <- weights*m*res
    U[is.na(U)] <- 0 
    
    #---bread---
    
    #summary(fit)$cov.unscaled is weighted
    I <- -solve(summary(fit)$cov.unscaled)/n  
    
  }
  if(inherits(x=fit, what="ah")){
  
    #---meat---
    
    #residuals are weighted 
    res <- predict(object=fit, type="residuals")
    rownames(res) <- fit$incl
    colnames(res) <- names(fit$coefficients)
    res <- expand(res, rownames(data))
    U <- res
    
  }
  if(inherits(x=fit, what="coxph")){
  
    #---meat for regression coefficients---
    
    #score residuals are unweighted, but computed with estimated coefficients 
    #from weighted model
    res <- residuals(fit, type="score")
    res <- expand(res, rownames(data)) 
    Ucoef <- weights*res   
    
    #---bread for regression coefficients---
    
    #vcov(fit) is weighted 
    Icoef <- -solve(vcov(fit))/n 
    
    if(missing(t)){
    
      U <- Ucoef
      colnames(U) <- names(fit$coef)
      U[is.na(U)] <- 0 
      I <- Icoef
    
    }else{
    
      nt <- length(t)
    
      #---meat and bread for baseline hazard---
    
      #check if left truncation
      varsLHS <- all.vars(fit$formula[[2]]) 
      t2 <- varsLHS[length(varsLHS)-1]
      if(length(varsLHS)==3)
        t1 <- varsLHS[1]  
      else
        t1 <- NULL
      time <- fit.detail$time
      #nevent is unweighted
      nevent <- fit.detail$nevent
      #hazard is weighted
      dH <- fit.detail$hazard
      names(dH) <- time
      #varhaz is dH/mean(exp(x*b)) where mean is in risk set
      #varhaz is weighted
      dHvar <- fit.detail$varhaz
      Hvar <- stepfun(time, c(0, cumsum(dHvar))) 
      p <- predict(object=fit, type="risk")
      m <- model.matrix(fit)
      names(p) <- rownames(m)
      p <- expand(p, rownames(data))  
      ncoef <- length(fit$coef) 
      #fit.detail$means is weighted, but not relative to the mean covariate 
      #in the sample, like all other outputs from coxph.detail,
      #subtracting fit$means fixes this
      means <- as.matrix(fit.detail$means)
      means <- means-matrix(fit$means, nrow=nrow(means), ncol=ncol(means), 
        byrow=TRUE)
      means <- expand(means, data[, t2])
      UH <- matrix(nrow=n, ncol=nt)
      IH <- matrix(nrow=nt, ncol=ncoef)
      for(j in 1:nt){
        #dividing with nevent accounts for ties,
        #but this assumes that H is computed with Breslow method for ties,
        #and that weights are equal within ties. 
        tmp1 <- n*expand(dH/nevent*(time<=t[j]), data[, t2]) 
        tmp1[is.na(tmp1)] <- 0
        tmp2 <- n*Hvar(pmin(t[j], data[, t2]))*p
        if(!is.null(t1))
          tmp2 <- tmp2-n*Hvar(data[, t1])*(data[, t1]<t[j])*p 
        UH[, j] <- tmp1-tmp2 
        dH.dbeta <- means*tmp1
        dH.dbeta[is.na(dH.dbeta)] <- 0
        IH[j, ] <- -colMeans(dH.dbeta)  
      }
      U <- cbind(Ucoef, UH) 
      colnames(U) <- c(names(fit$coef), paste0("H", t))
      U[is.na(U)] <- 0
      I <- rbind(cbind(Icoef, matrix(0, nrow=ncoef, ncol=length(t))), 
          cbind(IH, -diag(length(t))))
    }
      
  }
  if(inherits(x=fit, what="survfit")){ 
  
    #---meat---
    
    #check if left truncation
    #survfit object has no formula element, so need to get it from call,
    #need to use eval, since the fit$call$formula will be literary what the user
    #gave as argument, e.g. if formula=f, then fit$call$formula is f, not the 
    #formula contained in f  
    varsLHS <- all.vars(eval(fit$call$formula)[[2]]) 
    t2 <- varsLHS[length(varsLHS)-1]
    if(length(varsLHS)==3)
      t1 <- varsLHS[1]  
    else
      t1 <- NULL   
    #need to use summary(fit), since n.events and n.risk from fit 
    #behave strange when there is left truncation     
    ss <- summary(fit)
    strata <- ss$strata
    #n.strata is unweighted
    n.strata <- summary(strata)
    K <- length(n.strata)
    names.strata <- names(n.strata)
    time <- ss$time
    #n.event and n.risk are weighted
    n.event <- ss$n.event
    n.risk <- ss$n.risk
    dH <- n.event/n.risk
    names(dH) <- paste(time, strata)
    dHvar <- dH/n.risk
    #survfit object has no formula element, so need to get it from call,
    #need to use eval, since the fit$call$formula will be literary what the user
    #gave as argument, e.g. if formula=f, then fit$call$formula is f, not the 
    #formula contained in f
    vars <- all.vars(eval(fit$call$formula)[[3]])
    #note: strata is a function in the survival package
    strata.all <- strata(data[, vars, drop=FALSE]) 
    tmp1 <- matrix(nrow=n, ncol=K) 
    U <- matrix(nrow=n, ncol=K) 
    colnames(U) <- names.strata
    breaks <- c(0, cumsum(n.strata)) 
    for(k in 1:K){
      incl <- (breaks[k]+1):breaks[k+1] 
      Hvar <- stepfun(time[incl], c(0, cumsum(dHvar[incl]))) 
      #dividing with nevent[incl] account for ties, 
      #but this assumes that H is computed with Breslow method for ties,
      #and that weights are equal within ties.
      #multiplying with weights corrects for nevent being weighted; 
      #here we just want to divide with the actual number of events to account
      #for ties, not the weighted number of events
      tmp1.time <- n*dH[incl]/n.event[incl]*(time[incl]<=t)
      tmp1[, k] <- tmp1.time[match(paste(data[, t2], strata.all), 
        names(tmp1.time))]*weights
      tmp1[is.na(tmp1[, k]), k] <- 0
      sk <- names.strata[k]
      incl <- which(strata.all==sk)  
      tmp2 <- n*Hvar(pmin(t, data[incl, t2]))
      if(!is.null(t1))
        tmp2 <- tmp2-n*Hvar(data[incl, t1])*(data[incl, t1]<t)
      U[incl, k] <- tmp1[incl, k]-tmp2
    }
    
    #---bread---
    
    I <- diag(-1, K) 
    rownames(I) <- names.strata
    colnames(I) <- names.strata
    
  }
  
  U[is.na(U)] <- 0 
  return(list(I=I, U=U))

}

scs <- function(X, time, status, data, fitZ.L=NULL, max.time, max.time.beta, 
  n.sim){
  
  #stuff related to the model for Z
  Z <- all.vars(fitZ.L$formula)[1]
  IZ <- summary(object=fitZ.L)$cov.unscaled
  resZ <- residuals(object=fitZ.L, type="response")
  designZ <- model.matrix(object=fitZ.L)
  eps.theta <- (designZ*resZ)%*%IZ  
  g <- family(fitZ.L)$mu.eta
  dmu.deta <- g(predict(object=fitZ.L))  
  deta.dbeta <- designZ  
  dmu.dbeta <- dmu.deta*deta.dbeta
  E.dot <- dmu.dbeta
  p.dim <- dim(eps.theta)[2]
  Z.c <- data[, Z]-fitted(object=fitZ.L)
  
  #estimation of the target parameter
  stime1 <- sort(data[, time])
  n <- length(stime1)
  status_new <- data[, status]
  status_new[data[, time]>max.time] <- 0
  stime <- sort(data[status_new==1, time])
  res <- list()
  
  arvid1 <- data[, time] 
  arvid2 <- data[, X]
  res.tmp <- B_est1(arvid1, status_new, stime, stime1, Z.c, arvid2,
    max.time.beta, eps.theta, E.dot, p.dim)
  eps <- res.tmp$eps_B
  k <- length(stime<max.time.beta)
  deps <- eps
 
  deps[2:k,] <- eps[2:k, ]-eps[1:(k-1), ]
  eps.beta <- colSums(deps*res.tmp$at_risk)/res.tmp$tot_risk
  se_beta <- sqrt(sum(eps.beta^2))
  
  chisq_beta <- (res.tmp$beta/se_beta)^2
  pval_beta <- 1-pchisq(chisq_beta, df=1)
  ant.resamp <- n.sim
  GOF.resam0 <- matrix(0, ant.resamp, k)
  max.proc0 <- numeric(ant.resamp)
  GOF.resam <- matrix(0, ant.resamp, k)
  max.proc <- CvM.proc <- numeric(ant.resamp) 
  eps.const.eff <- eps[1:k, ]

  for (j in 1:k)
    eps.const.eff[j, ] <- eps[j,]-stime[j]*eps.beta
  for (j1 in 1:ant.resamp){
    G <- rnorm(n, 0, 1)
    tmp.mat0 <- eps%*%matrix(G, n, 1)
    GOF.resam0[j1,] <- c(tmp.mat0)
    max.proc0[j1] <- max(abs(GOF.resam0[j1, ]))
    tmp.mat <- eps.const.eff%*%matrix(G, n, 1)
    GOF.resam[j1,] <- c(tmp.mat)
    max.proc[j1] <- max(abs(GOF.resam[j1, ]))
    CvM.proc[j1] <- sum(GOF.resam[j1, ]^2*c(stime[2:k]-stime[1:(k-1)], 
      max.time.beta-stime[k]))
  }
  GOF.proc0 <- res.tmp$B[1:k]
  max.obs0 <- max(abs(GOF.proc0))
  GOF.proc <- res.tmp$B[1:k]-res.tmp$beta*stime
  max.obs <- max(abs(GOF.proc))
  CvM.obs <- sum(GOF.proc^2*c(stime[2:k]-stime[1:(k-1)],max.time.beta-stime[k]))
  pval_0 <- sum(max.proc0>max.obs0)/ant.resamp
  pval.sup <- sum(max.proc>max.obs)/ant.resamp
  pval.CvM <- sum(CvM.proc>CvM.obs)/ant.resamp
  res$stime <- res.tmp$stime
  res$B <- res.tmp$B
  res$se_B <- res.tmp$se
  res$pval_0 <- pval_0
  res$eps_B <- res.tmp$eps_B
  est <- res.tmp$beta
  names(est) <- X
  res$est <- est
  vcov <- matrix(se_beta^2)
  rownames(vcov) <- X
  colnames(vcov) <- X
  res$vcov <- vcov
  res$pval_psi <- pval_beta
  res$pval_GOF_sup <- pval.sup
  res$pval_GOF_CvM <- pval.CvM
  res$fitZ.L <- fitZ.L
  if (n.sim>50)
    res$GOF.resamp <- rbind(stime[1:k], GOF.proc, GOF.resam[1:50,])
  return(res)

}

summary.ivbounds <- function(object, ...) {
  p0min <- object$p0["min"]
  p0max <- object$p0["max"]
  p1min <- object$p1["min"]
  p1max <- object$p1["max"]
  CRDmin <- p1min-p0max
  CRDmax <- p1max-p0min
  CRRmin <- p1min/p0max
  CRRmax <- p1max/p0min
  CORmin <- p1min/(1-p1min)/(p0max/(1-p0max))
  CORmax <- p1max/(1-p1max)/(p0min/(1-p0min))
  lwr <- c(p0min, p1min, CRDmin, CRRmin, CORmin)
  upr <- c(p0max, p1max, CRDmax, CRRmax, CORmax)
  coef.table <- cbind(lwr, upr)
  rownames(coef.table) <- c("p0", "p1", "CRD", "CRR", "COR")
  colnames(coef.table) <- c("lower", "upper")
  p0.symbolicmin <- object$p0.symbolic["min"]
  p0.symbolicmax <- object$p0.symbolic["max"]
  p1.symbolicmin <- object$p1.symbolic["min"]
  p1.symbolicmax <- object$p1.symbolic["max"]
  lwr.symbolic <- c(p0.symbolicmin, p1.symbolicmin)
  upr.symbolic <- c(p0.symbolicmax, p1.symbolicmax)
  coef.table.symbolic <- cbind(lwr.symbolic, upr.symbolic)
  rownames(coef.table.symbolic) <- c("p0", "p1")
  colnames(coef.table.symbolic) <- c("lower", "upper")
  ans <- list(call=object$call, coefficients=coef.table,
    coefficients.symbolic=coef.table.symbolic, 
    IVinequality=object$IVinequality, conditions=object$conditions)
  class(ans) <- "summary.ivbounds"
  return(ans)
}


summary.ivmod <- function(object, ...) {
  s.err <- sqrt(diag(as.matrix(object$vcov)))
  zvalue <- object$est/s.err
  pvalue <- 2*pnorm(-abs(zvalue))
  coef.table <- as.matrix(cbind(object$est, s.err, zvalue, 
        pvalue))
  dimnames(coef.table) <- list(names(object$est), c("Estimate", 
      "Std. Error", "z value", "Pr(>|z|)"))
  ans <- list(call=object$call, coefficients=coef.table, t=object$t)
  if(inherits(x=object, what="ivah") & object$input$estmethod=="g"){
    test0 <- cbind(object$pval_0) 
    colnames(test0)  <- c("Supremum-test pval")
    rownames(test0)<- c(object$input$X) 
    test_gof <- cbind(object$pval_GOF_sup)
    colnames(test_gof) <- c("Supremum-test pval")
    rownames(test_gof) <- c("        ")
    ans$test0 <- test0
    ans$test_gof <- test_gof
  }
  class(ans) <- "summary.ivmod"
  return(ans)
}

tsest <- function(fitX.LZ, fitY.LX, data, ctrl=FALSE, clusterid=NULL, 
  vcov.fit=vcov.fit){
  
  #preliminaries
  n <- nrow(data)
  weights <- expand(fitX.LZ$prior.weights, rownames(data))
  if(!is.null(clusterid))
    ncluster <- length(unique(data[, clusterid]))
  
  #collect features of fitX.LZ
  X <- as.character(fitX.LZ$formula[2])
  resX <- expand(residuals(object=fitX.LZ, type="response"), rownames(data))
  designX <- expand(model.matrix(object=fitX.LZ), rownames(data))
  nX <- length(fitX.LZ$coef)
   
  #estimate
  Xhat <- expand(predict(object=fitX.LZ, type="response"), rownames(data))
  data.Xhat <- data
  data.Xhat[, X] <- Xhat
  call <- fitY.LX$call
  frml <- fitY.LX$formula
  if(ctrl){
    data.Xhat[, "R"] <- resX
    call$formula <- update(frml, .~.+R)
  }
  call$data <- data.Xhat 
  fitY.LX <- eval(call, envir=environment(frml))
  est <- fitY.LX$coef 
  
  #collect features of fitY.LX
  nY <- length(fitY.LX$coef)  
  
  #converged?
  if(inherits(x=fitY.LX, what="glm")){  
    converged <- fitY.LX$converged
  }                               
  else{                           
    converged <- TRUE             
  }     
   
  #compute variance if requested, and if solution to estimating equation was found
  if(vcov.fit & converged){
    sandwich.fitX.LZ <- sandwich(fit=fitX.LZ, data=data, weights=weights)
    sandwich.fitY.LX <- sandwich(fit=fitY.LX, data=data, weights=weights)
    UX <- sandwich.fitX.LZ$U
    UY <- sandwich.fitY.LX$U
    U <- cbind(UX, UY)
    if(!is.null(clusterid))
      U <- aggr(x=U, clusters=data[, clusterid])
    J <- var(U)  
    IX <- cbind(sandwich.fitX.LZ$I, matrix(0, nrow=nX, ncol=nY)) 
    UYfunc <- function(b){
      Xhat <- as.vector(family(fitX.LZ)$linkinv(designX%*%b)) 
      data.Xhat[, X]  <- Xhat
      if(ctrl)
        data.Xhat[, "R"] <- data[, X]-Xhat
      if(inherits(x=fitY.LX, what="glm")){
        designY <- expand(model.matrix(object=fitY.LX$formula, data=data.Xhat), 
          rownames(data))
        Yhat <- as.vector(family(fitY.LX)$linkinv(designY%*%est))
        Y <- expand(fitY.LX$y, rownames(data))
        resY <- Y-Yhat
        UY <- weights*resY*designY
      }
      if(inherits(x=fitY.LX, what="coxph")){
        fitY.LX$call$data <- data.Xhat
        resY <- expand(residuals(object=fitY.LX, type="score"), rownames(data))
        UY <- as.matrix(weights*resY)
      }
      if(inherits(x=fitY.LX, what="ah")){
        fitY.LX$data$X <- model.matrix(object=fitY.LX$formula, 
          data=data.Xhat)[, -1, drop=FALSE] 
        resY <- predict(object=fitY.LX, type="residuals")
        rownames(resY) <- fitY.LX$incl
        colnames(resY) <- names(fitY.LX$coefficients)
        resY <- expand(resY, rownames(data))
        UY <- resY
      }
      UY[is.na(UY)] <- 0  
      colMeans(UY)
    } 
    if(inherits(x=fitY.LX, what="glm"))
      IY <- cbind(jacobian(func=UYfunc, x=fitX.LZ$coef),
        -solve(summary(object=fitY.LX)$cov.unscaled)/n)
    if(inherits(x=fitY.LX, what="coxph"))
      IY <- cbind(jacobian(func=UYfunc, x=fitX.LZ$coef),
        -solve(vcov(object=fitY.LX))/n)   
    if(inherits(x=fitY.LX, what="ah"))
      IY <- cbind(jacobian(func=UYfunc, x=fitX.LZ$coef),
        -fitY.LX$D/n)
    I <- rbind(IX, IY)
    rownames(I) <- colnames(U)
    colnames(I) <- colnames(U)
    if(is.null(clusterid))
      vcov <- (solve(I)%*%J%*%t(solve(I))/n)
    else
      vcov <- (solve(I)%*%J%*%t(solve(I))*ncluster/n^2)
    vcov <- vcov[(nX+1):(nX+nY), (nX+1):(nX+nY), drop=FALSE]
  }
  else{
    U <- NA
    I <- NA
    vcov <- matrix(NA, nrow=nY, ncol=nY)  
  }
  rownames(vcov) <- names(est)
  colnames(vcov) <- names(est)
                              
  result <- list(est=est, vcov=vcov, estfunall=U, d.estfun=I, fitY.LX=fitY.LX, 
    converged=converged)
  class(result) <- "tsest" 

  return(result)
  
}

#computes estimating function for each subject
Upsifun <- function(b, Z, X, Y, type, fitZ.L, fitX.LZ, fitX.L, fitY.LZX, 
  npsi, nZ, nX.LZ, nX.L, nY, designpsi, designZ, designX.LZ, designX.L, designY, 
  weights, data){

  n <- nrow(data)
  psi <- b[1:npsi]
  g <- designpsi*data[, X]
  lppsi <- as.vector(g%*%psi)
  if(is.null(fitX.LZ)){
    nd <- nZ
  }else{
    nd <- nX.LZ+nX.L  
  }
  if(type=="identity")
    Y0hat <- data[, Y]-lppsi
  if(type=="log")
    Y0hat <- data[, Y]*exp(-lppsi)
  if(type=="logit"){
    bY <- b[(npsi+nd+1):(npsi+nd+nY)]
    lpY <- as.vector(designY%*%bY)
    Y0hat <- plogis(lpY-lppsi)  
  }
  if(type=="coxph"){
    if(inherits(x=fitY.LZX, what="survfit")){
        H <- b[(npsi+nd+1):(npsi+nd+nY)]
        #survfit object has no formula element, so need to get it from call,
        #need to use eval, since the fit$call$formula will be literary what the user
        #gave as argument, e.g. if formula=f, then fit$call$formula is f, not the 
        #formula contained in f
        vars <- all.vars(eval(fitY.LZX$call$formula)[[3]])
        strata.all <- strata(data[, vars, drop=FALSE])
        Y0hat <- exp(-H[strata.all]*exp(-lppsi))  
      }
    if(inherits(x=fitY.LZX, what="coxph")){
      bY <- b[(npsi+nd+1):(npsi+nd+nY-1)]
      H0 <- b[npsi+nd+nY]
      lpY <- as.vector(designY%*%bY)  
      Y0hat <- exp(-H0*exp(lpY-lppsi))
    }
  }
  if(is.null(fitX.LZ)){
    bZ <- b[(npsi+1):(npsi+nZ)]
    d <- designpsi*(data[, Z]-as.vector(family(fitZ.L)$linkinv(designZ%*%bZ)))     
  }else{
    bX.LZ <- b[(npsi+1):(npsi+nX.LZ)]
    bX.L <- b[(npsi+nX.LZ+1):(npsi+nX.LZ+nX.L)] 
    d <- designpsi*(as.vector(family(fitX.LZ)$linkinv(designX.LZ%*%bX.LZ))-
      as.vector(family(fitX.L)$linkinv(designX.L%*%bX.L)))      
  } 
  Upsi <- weights*d*Y0hat
  Upsi[is.na(Upsi)] <- 0
  colnames(Upsi) <- names(psi)

  return(Upsi)
  
}

#extracts and computes various useful stuff for G-estimation  
utilityfun <- function(Z, X, Y, type, data, formula, y, fitY.LZX, fitX.LZ, 
  fitX.L, fitZ.L, fit.detail){
 
  n <- nrow(data)
  
  #stuff related to Z or X
  if(is.null(fitX.LZ)){
    nX.LZ <- NULL
    designX.LZ <- NULL
    nX.L <- NULL
    designX.L <- NULL
    nZ <- length(fitZ.L$coef)
    designZ <- expand(model.matrix(object=fitZ.L), rownames(data)) 
    weights <- expand(fitZ.L$prior.weights, rownames(data))
    est.nuisance <- fitZ.L$coef 
  }else{ 
    nZ <- NULL
    designZ <- NULL
    nX.LZ <- length(fitX.LZ$coef)
    designX.LZ <- expand(model.matrix(object=fitX.LZ), rownames(data)) 
    weights <- expand(fitX.LZ$prior.weights, rownames(data))
    nX.L <- length(fitX.L$coef)
    designX.L <- expand(model.matrix(object=fitX.L), rownames(data))
    est.nuisance <- c(fitX.LZ$coef, fitX.L$coef)  
  }
  
  #stuff related to Y 
  t1 <- NULL  
  if(type=="identity" | type=="log"){
    nY <- 0  
    designY <- NULL
  }
  if(type=="logit"){
    nY <- length(fitY.LZX$coef) 
    designY <- expand(model.matrix(object=fitY.LZX), rownames(data))
    est.nuisance <- c(est.nuisance, fitY.LZX$coef)
  }
  if(type=="coxph"){ 
    if(inherits(x=fitY.LZX, what="survfit")){
      #survfit object has no formula element, so need to get it from call,
      #need to use eval, since the fit$call$formula will be literary what the user
      #gave as argument, e.g. if formula=f, then fit$call$formula is f, not the 
      #formula contained in f
      temp <- all.vars(eval(fitY.LZX$call$formula)[[2]])
      designY <- NULL
      H <- Hfun(fit=fitY.LZX, data=data)
      nY <- length(H)
      Hy <- vector(length=nY)
      for(k in 1:nY)
        Hy[k] <- H[[k]](y)
      est.nuisance <- c(est.nuisance, Hy)
    }
    if(inherits(x=fitY.LZX, what="coxph")){
      temp <- all.vars(fitY.LZX$formula[[2]])
      designY <- expand(model.matrix(object=fitY.LZX), rownames(data))
      H <- Hfun(fit=fitY.LZX, data=data, fit.detail=fit.detail)  
      est.nuisance <- c(est.nuisance, fitY.LZX$coef, H(y))  
      nY <- length(fitY.LZX$coef)+1      
    } 
    #left-truncation?    
    if(length(temp)==3) 
      t1 <- temp[1]    
  }
  
  #stuff related to psi 
  temp <- model.matrix(object=formula, data=data)
  #If formula==~1, then model.matrix does not use the row names in data!
  #Note: cannot remove if(formula==~1), since rownames(temp) <- rownames(data)
  #does not work if missing, since then nrow(temp)!=nrow(data). However,
  #if formula==~1, then missing will not reduce size of temp. 
  if(formula==~1)
    rownames(temp) <- rownames(data)
  designpsi <- expand(temp, rownames(data))
  npsi <- ncol(designpsi)   
  
  return(list(est.nuisance=est.nuisance, npsi=npsi, nZ=nZ, nX.LZ=nX.LZ,
    nX.L=nX.L, nY=nY, designpsi=designpsi, designZ=designZ, 
    designX.LZ=designX.LZ, designX.L=designX.L, designY=designY, 
    weights=weights, t1=t1))

}






