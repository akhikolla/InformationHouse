augSIMEX<-function(mainformula = formula(data), mismodel = pi|qi~1, meformula = NULL, family = gaussian,
                   data,validationdata,
                   err.var, mis.var, err.true, mis.true, err.mat = NULL, cppmethod = TRUE,
                   repeated = FALSE, repind=list(),
                   subset,  offset, weights, na.action, scorefunction=NULL,
                   lambda = NULL, M = 5, B = 20, nBoot = 50, extrapolation = c("quadratic","linear"),
                   bound = 8, initial = NULL,...)

{ call <- match.call()

  if (!missing(family)){
    ### makesure family is matched
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    
    ### Obtain the function in the family object
    variance <- family$variance
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
      stop("'family' argument seems not to be a valid family object",
           call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    
    unless.null <- function(x, if.null) if (is.null(x))
      if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    
    if (family$family=="gaussian") {dispersionpar = TRUE} else {
      dispersionpar = FALSE
    }
    ### check if fast approach available
    cppmethod<-fastapproach(family,cppmethod)$fast
    scorefun<-fastapproach(family,cppmethod)$scorefun
  } else {
    if (is.null(scorefunction)) stop("Family is missing and scorefunction is not specified. ")
    sfun <-scorefunction
    scorefun<-score.modifieduser
  }

 

  ### setup the weights and offset
  if (missing(weights)) weights<-NULL
  if (missing(offset))  offset<-NULL

  if ((repeated) & (length(repind)==0) ){
     stop("The repind should be specified for each repeated covariates.")}

  ### make data into design matrix
  if (missing(data)) data <- environment(mainformula)
  if (repeated) {
    if (!missing(repind)) {
      for (i in 1:length(names(repind))){
        if (!names(repind)[i] %in% names(data)) {
          data<-cbind(data,rep(0,dim(data)[1]))
          colnames(data)[dim(data)[2]]<-names(repind)[i]}
        }
    }
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("mainformula", "data", "subset","weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  if (!repeated) {names(mf)[2]<-"formula"} else {
    names(mf)[2]<-"formula"
    mf$data <- data
  }
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  xlevels <- .getXlevels(mt, mf)
  Response <- model.response(mf, "any")
  nsize<-length(Response)
  if (length(dim(nsize)) == 1L) {
    nm <- rownames(nsize)
    dim(nsize) <- NULL
    if (!is.null(nm))
      names(nsize) <- nm
  }
  mainCovariates <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts) else matrix(, nsize, 0L)
  contrasts <- attr(mainCovariates, "contrasts")
  if (all(mainCovariates[,1]==1)) intercept<-TRUE
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != nsize)
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), nsize), domain = NA)
  }else {offset<- rep.int(0, nsize)}
  if (is.null(weights)) weights <- rep.int(1, nsize)
   else { if (length(weights) != nsize)
     stop(gettextf("number of weights is %d should equal %d (number of observations)",
                       length(weights), nsize), domain = NA)}
  qr <- qr(mainCovariates)

  ### rearrange the validationData
  val.covariates<-c(mis.var,err.var,err.true, mis.true)
  if (missing(validationdata)){ validationdata<-mget(unique(val.covariates),envir=parent.frame())}
  if (!repeated){
    if (!all(val.covariates %in% names(validationdata)))
      {stop(paste("Variable",
                   val.covariates[!val.covariates %in% names(validationdata)],
                   "is not included in the validation data."))}}
   else {if (!all(c(mis.var,mis.true) %in% names(validationdata))) {
     stop(paste("Variable",
                mis.var, "or", mis.true,
                "is not included in the validation data."))
   }}
  colname<-colnames(data)

  ## SIMEXvariable
  err.var = unique(err.var)
  nerr.var = length(err.var)
  mis.var = unique(mis.var)
  nmis.var = length(mis.var)
  if (!is.character(err.var) | nerr.var > length(colname)) {
    stop("Invalid err.var object")
  }
  if (!is.character(mis.var) | nmis.var > length(colname)) {
    stop("Invalid mis.var object")
  }
  if (!all(err.var %in% colname)) {
    stop("Error-prone variable must be selected from the data")
  }
  if (!all(mis.var %in% colname)) {
    stop("Misspecified variable must be included in the data")
  }
  if (!(repeated == FALSE) & !(repeated == TRUE)) {
    stop("Repeated indicator should only be 'TRUE' or 'FALSE'. ")
  }
  if (!is.null(err.mat)) {
    err.mat = as.matrix(err.mat)
    if (!is.numeric(err.mat) | any(err.mat < 0)) {
      stop("Invalide err.mat object, err.mat must be a square symmetric numeric matrix")
    }
    if (nrow(err.mat) != ncol(err.mat)) {
      stop("err.mat must be a square matrix")
    }
    if (length(err.var) != nrow(err.mat)) {
      stop("SIMEXvariable and err.mat have non-conforming size")
    }
    SSigma <- err.mat
    dimnames(SSigma) <- NULL
    if (!isTRUE(all.equal(SSigma, t(SSigma)))) {
      warning("err.mat is numerically not symmetric")
    }
  }


  ###SIMEX parameters
  if (length(B) != 1) {
    stop("B must be positive integer")
  }
  if (!is.numeric(B) | B <= 0) {
    stop("B must be positive integer")
  }
  else {
    B = ceiling(B)
  }
  if (is.null(lambda)) {lambda<-seq(from=0,to=2,length.out=M)}
  if (!is.vector(lambda) | !is.numeric(lambda)) {
    stop(":Invalide lambda object")
  }
  if (any(lambda < 0)) {
    warning("Lambda should be positive values. Negative values will be ignored",
            call. = FALSE)
    lambda <- lambda[lambda >= 0]
  }

  # extrapolation<-match.arg(extrapolation)

  ### Step 1: Simulation step
  temp.Results1<-Getalpha(validationdata,as.Formula(mismodel),mis.var,mis.true,err.var)
  alphahat1<-temp.Results1$alphahat1
  alphahat2<-temp.Results1$alphahat2

  ### estimate the empirical covariance matrix
  if (!repeated){  
    if (missing(err.mat)|is.null(err.mat)){
      if (length(meformula)==0) {
      Models_res<-lapply(1:length(err.var),FUN = function(i){
        model<-lm(validationdata[,err.var[i]]~-1+offset(validationdata[,err.true[i]]),data=data.frame(validationdata))
        return(model$residuals)
      })} else {
        meformula = as.Formula(meformula)
        if (length(meformula)[1]!=length(err.var)) {stop(gettextf("the length of meformula responses (left hand side) is %d and should be equal to the number of error-prone covariates %d",
                                                  length(meformula), length(err.var)), domain = NA)}
      resname.me <- all.vars(meformula)[1:length(meformula)[1]]
      me.mf <- model.frame(meformula, data = validationdata, na.action = na.omit)
      
      index<-match(err.var,resname.me)
      if (length(index)!=length(err.var)) {stop("incorrect measurement error model specified")}
      if (length(meformula)[2]==1) {
          Models_res<-lapply(index,FUN = function(i){
            model<-lm(formula(meformula,lhs=i,rhs=1),data = me.mf)
            return(model$residuals)
          })
        } else {
          Models_res<-lapply(index,FUN = function(i){
            model<-lm(formula(meformula,lhs=i,rhs=i),data = me.mf)
            return(model$residuals)
          })
        }
      }
      Models_res<-matrix(unlist(Models_res),ncol=length(err.var))
      Sigma_ehat<-cov(Models_res)
    } else Sigma_ehat<-cov(err.mat)
  }
  
  ### change the Xstar in mainCovariates into X, Zstar into Z
  changeindex<-match(err.var,colnames(mainCovariates))
  colnames(mainCovariates)[changeindex]<-err.true
  colnames(mainCovariates)[which(colnames(mainCovariates)==mis.var)]<-mis.true

  main.new<-mainCovariates
  Wnames<- setdiff(colnames(main.new),c(val.covariates,unlist(repind)))
  Wmatrix<-main.new[,Wnames]
  NSam<-dim(main.new)[1]
  nbeta<-dim(mainCovariates)[2]
  
  if (dispersionpar) {
    nbeta <- nbeta + 1
    if (is.null(initial)) { initialv <- c(rep(0,nbeta-1),0.25)} else{
      if (length(initial)!=nbeta) {stop(paste0("The length of initial values should be ",nbeta,"."))}
      initialv <- initial 
     }
  } else {
    if (is.null(initial)) {initialv <- rep(0,nbeta)} else {
      if (length(initial)!=nbeta) {stop(paste0("The length of initial values should be ",nbeta,"."))}
      initialv <- initial 
    }
    }
  
  if (repeated) {
      if (options("na.action")=="na.omit") {
        completeindex <- complete.cases(data[,! colname %in% unlist(repind)])
        imputeData <- data[completeindex,] } else imputeData<-data 
  } else imputeData<-mainCovariates

  betahat<-lapply(1:M,FUN=function(i){
    betab_lam<-lapply(1:B,FUN=function(x){
      X.impute<-imputeX(imputeData,err.true,err.var,lambda[i],Sigma_ehat,nsize,repeated,repind)
      main.new[,err.true]<-X.impute
      main.new.df<-data.frame(main.new)

      phat0<-predict(alphahat1,newdata=main.new.df,type="response")
      qhat0<-predict(alphahat2,newdata=main.new.df,type="response")

      DataM <- cbind(X.impute,Wmatrix,main.new[,mis.true])
      DataM0 <- cbind(X.impute,Wmatrix,rep(0,nsize))
      DataM1 <- cbind(X.impute,Wmatrix,rep(1,nsize))

      if (is.null(scorefunction)){
          if (cppmethod) {betasolve<-nleqslv(initialv,scorefun,Y=Response,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights,offset=offset)}
          else {betasolve<-nleqslv(initialv,scorefun,Y=Response,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights,offset=offset,linkinv=linkinv,var=variance,mueta=mu.eta)}
      } else {betasolve<-nleqslv(initialv,scorefun,Y=Response,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights,offset=offset,sfun=sfun)}
      betareturn <- betasolve$x
      
      

      if (any(abs(betareturn)>bound, na.rm = T)) { return(rep(NA,nbeta*2))}

      return(betareturn)
    })
    
    betahat_lam_M<-matrix(unlist(betab_lam),ncol=nbeta,byrow=T)
    betahat_lam<-colMeans(betahat_lam_M[,1:nbeta],na.rm=T)
    return(c(betahat_lam))
  })

  betamatrix<-matrix(unlist(betahat),ncol=nbeta,byrow=T)
  coefs <- extrapolate(extrapolation,betamatrix,lambda,nbeta)
  if (dispersionpar){
    if (extrapolation=="both"){coefs <- coefs[c(1:(nbeta-1),(nbeta+1):(nbeta*2-1))]}
    else {coefs <- coefs[c(1:(nbeta-1))]}
    betamatrix <- betamatrix[,c(1:(nbeta-1))]
  }
  
  ### bootstrap procedure
  betahat_boot_all<-lapply(1:nBoot,FUN=function(t){

    sample.boot<-sample(NSam,replace=T)
    imputeData_boot<-imputeData[sample.boot,]
    main.new.boot<-main.new[sample.boot,]
    Wmatrix_b<-main.new.boot[,Wnames]
    Response_b<-Response[sample.boot]
    
    betahat_boot<-lapply(1:M,FUN=function(i){
      betab_lam<-lapply(1:B,FUN=function(x){
        X.impute<-imputeX(imputeData_boot,err.true,err.var,lambda[i],Sigma_ehat,nsize,repeated,repind)
        options(warn=-1)
        main.new.boot[,err.true]<-X.impute
        main.new.df<-data.frame(main.new.boot)
        
        phat0<-predict(alphahat1,newdata=main.new.df,type="response")
        qhat0<-predict(alphahat2,newdata=main.new.df,type="response")
        
        DataM <- cbind(X.impute,Wmatrix_b,main.new.boot[,mis.true])
        DataM0 <- cbind(X.impute,Wmatrix_b,rep(0,nsize))
        DataM1 <- cbind(X.impute,Wmatrix_b,rep(1,nsize))
        options(warn=0)
        
        
        if (is.null(scorefunction)){
            if (cppmethod) {betasolve<-nleqslv(initialv,scorefun,Y=Response_b,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights[sample.boot],offset=offset[sample.boot])}
            else {betasolve<-nleqslv(initialv,scorefun,Y=Response_b,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights[sample.boot],offset=offset[sample.boot],linkinv=linkinv,var=variance,mueta=mu.eta)}
        } else {betasolve<-nleqslv(initialv,scorefun,Y=Response_b,DataM=DataM,DataM0=DataM0,DataM1=DataM1,phat=phat0,qhat=qhat0,weight=weights[sample.boot],offset=offset[sample.boot],sfun=sfun)}
        betareturn <- betasolve$x
        
        if (any(abs(betareturn)>bound,na.rm=T)) { return(rep(NA,nbeta*2))}
        
        return(betareturn)
      })
      
      betahat_lam_M<-matrix(unlist(betab_lam),ncol=nbeta,byrow=T)
      betahat_lam<-colMeans(betahat_lam_M[,1:nbeta],na.rm=T)
      return(c(betahat_lam))
    })

    betamatrix_boot<-matrix(unlist(betahat_boot),ncol=nbeta,byrow=T)
    coefs <- tryCatch(extrapolate(extrapolation,betamatrix_boot,lambda,nbeta),error = function(e) {
      if (extrapolation=="both") return(rep(NA,2*nbeta)) else return(rep(NA,nbeta))})
     
    return(coefs)
  })
  
  if (extrapolation=="both"){
    if (dispersionpar){
      betahat_boot_all_M<-matrix(unlist(betahat_boot_all),
                                 ncol=nbeta*2,byrow=T)
      vcov<-apply(betahat_boot_all_M[,c(1:(nbeta-1),(nbeta+1):(nbeta*2-1))],MARGIN = 2,FUN=sd, na.rm=T)
      nbeta <- nbeta -1 
    } else{
      betahat_boot_all_M<-matrix(unlist(betahat_boot_all),ncol=nbeta*2,byrow=T)
      vcov<-apply(betahat_boot_all_M,MARGIN = 2,FUN=sd, na.rm=T)
    }
  } else{
    if (dispersionpar){
      betahat_boot_all_M<-matrix(unlist(betahat_boot_all),ncol=nbeta,byrow=T)
      vcov<-var(betahat_boot_all_M[,c(1:(nbeta-1))]) 
      nbeta <- nbeta -1 
    } else{
      betahat_boot_all_M<-matrix(unlist(betahat_boot_all),ncol=nbeta,byrow=T)
      vcov<-var(betahat_boot_all_M) 
    }
  }
  
  
  ### results wrap-up
  good <- weights > 0
  if (extrapolation=="both") {
    names(coefs)<-names(vcov)<-rep(c(err.true,Wnames,mis.true),2)
  } else {names(coefs)<-names(vcov)<-c(err.true,Wnames,mis.true)}

  if  (extrapolation=="both") { colnames(betamatrix)<-names(coefs)[1:nbeta]} else{
  colnames(betamatrix)<-c(names(coefs))}
  
  ### evaluate the performance from glm theory
  if (is.null(scorefunction)){
    eta <- drop(mainCovariates[,names(coefs)] %*% t(t(coefs)))
    mu <- linkinv(eta <- eta + offset)
    
    varmu <- variance(mu)[good]
    if (anyNA(varmu))
      stop("NAs in V(mu)")
    if (any(varmu == 0))
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good])))
      stop("NAs in d(mu)/d(eta)")
    
    residuals <- (Response - mu)/mu.eta(eta)
    dev<-sum(dev.resids(Response, mu, weights))
    wtdmu <- if (intercept) sum(weights * Response)/sum(weights) else linkinv(offset)
    nulldev <- sum(dev.resids(Response, wtdmu, weights))
    ## calculate df
    n.ok <- nsize - sum(weights==0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- length(coefs)
    resdf  <- n.ok - rank
    ## calculate AIC
    aic.model <- aic(Response, 0, mu, weights, dev) + 2*rank
  } else {
    eta <- NA
    mu <- NA
    varmu <- NA
    residuals <- NA
    dev<- NA
    wtdmu <- NA
    nulldev <- NA
    ## calculate df
    n.ok <- nsize - sum(weights==0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- length(coefs)
    resdf  <- n.ok - rank
    ## calculate AIC
    aic.model <- NA
  }
  
  
  ### Move the intercept in the front of parameters
  if ((intercept) && (extrapolation!="both")) {
    intindex <- which(names(coefs)=="(Intercept)")
    otherindex <- 1:nbeta
    otherindex <- otherindex[-intindex]
    coefs <- coefs[c(intindex,otherindex)]
    vcov <- vcov[c(intindex,otherindex),c(intindex,otherindex)]
    betamatrix <- betamatrix[,c(intindex,otherindex)]
  }


  output<-list(coefficients = coefs,
            vcov = vcov, fitted.values = mu,
            call = call, family = family,
            formula = mainformula, mismodel = mismodel, terms = mt,
            err.var = err.var, mis.var = mis.var, err.true = err.true, mis.true = mis.true, err.mat = err.mat,
            M = M, B = B, nBoot = nBoot, extrapolation = extrapolation, weights = weights,
            lambda = lambda, coefmatrix = betamatrix, linear.predictors = eta, deviance = dev, aic = aic.model,
            rank = rank, null.deviance = nulldev, df.residual = resdf, df.null = nulldf, residuals = residuals, 
            qr = qr, x = mainCovariates, y = Response, model = mf, contrasts = contrasts, xlevels = xlevels)

  class(output)<- "augSIMEX"

  return(output)

}


