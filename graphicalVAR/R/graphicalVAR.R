computePCC <- function(x)
{
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- as.matrix(Matrix::forceSymmetric(x))
  return(x)
}

computePDC <- function(beta,kappa){
  if (ncol(beta) == nrow(beta)+1){
    beta <- beta[,-1,drop=FALSE]
  }
  sigma <- solve(kappa)
  t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}

graphicalVAR <-
  function(
    data, # A n by p data frame containing repeated measures
    nLambda = 50, # Either single value or vector of two corresponding to c(kappa, beta)
    verbose = TRUE,
    gamma = 0.5,
    scale = TRUE,
    lambda_beta,
    lambda_kappa, maxit.in = 100, maxit.out = 100,
    deleteMissings = TRUE,
    penalize.diagonal = TRUE,
    lambda_min_kappa = 0.05,
    lambda_min_beta = lambda_min_kappa,
    mimic = c("current","0.1.2","0.1.4","0.1.5","0.2"),
    vars,
    beepvar,
    dayvar,
    idvar,
    lags = 1,
    centerWithin = TRUE,
    likelihood = c("unpenalized","penalized") # compute likelihood based on unpenalized (sparseTSCGM 2.3) or penalized (sparseTSCGM 2.5) contemporaneous effects
  ){
    
    mimic <- match.arg(mimic)
    if (mimic == "0.1.2"){
      if (lambda_min_beta != lambda_min_kappa){
        warning("mimic = 0.1.2 only uses lambda_min_kappa, not lambda_min_beta")
      }
      if (lambda_min_kappa != 0.01){
        warning("Set lambda_min_kappa = 0.01 to mimic 0.1.2 default behavior")
      }
    }
    
    # If list:
    if (is.list(data) && !is.data.frame(data)){
      if (!("data_c" %in% names(data) & "data_l" %in% names(data))){
        stop("'data_c' and 'data_l' must be contained in 'data'")
      }
      data_c <- data$data_c
      data_l <- data$data_l
    } else{
      if (mimic == "0.1.5"){
        # Check input:
        if (is.data.frame(data)){
          data <- as.matrix(data)
        }
        stopifnot(is.matrix(data))
        
        # Center data:
        data <- scale(data, TRUE, scale)
        
        # Compute current and lagged data:
        data_c <- data[-1,,drop=FALSE]
        data_l <- cbind(1,data[-nrow(data),,drop=FALSE])
      } else {
        data <- tsData(as.data.frame(data),
                       vars = vars,
                       beepvar = beepvar,
                       dayvar = dayvar,
                       idvar = idvar,
                       scale = scale,
                       centerWithin = centerWithin,
                       lags = lags)
        data_c <- data$data_c
        data_l <- data$data_l
      }
    }
    data_c <- as.matrix(data_c)
    data_l <- as.matrix(data_l)
    
    Nvar <- ncol(data_c)
    Ntime <- nrow(data_c)
    
    
    # Delete missing rows:
    if (any(is.na(data_c)) || any(is.na(data_l))){
      
      if (deleteMissings){
        
        warnings("Data with missings deleted")
        
        missing <- rowSums(is.na(data_c)) > 0 | rowSums(is.na(data_l)) > 0
        data_c <- data_c[!missing,]
        data_l <- data_l[!missing,]
        
      } else {
        stop("Missing data not supported")
      }
      
    }
    
    # Generate lambdas (from SparseTSCGM package):
    if (missing(lambda_beta) | missing(lambda_kappa)){
      if (mimic == "0.1.2"){
        lams <- SparseTSCGM_lambdas(data_l, data_c, nLambda, lambda.min.ratio=lambda_min_kappa)      
      } else {
        lams <- generate_lambdas(data_l, data_c, nLambda,nLambda, lambda_min_kappa=lambda_min_kappa,lambda_min_beta=lambda_min_beta,penalize.diagonal=penalize.diagonal,
                                 version0.1.4 = mimic == "0.1.4")      
      }
      if (missing(lambda_beta)){
        lambda_beta <- lams$lambda_beta
      }
      if (missing(lambda_kappa)){
        lambda_kappa <- lams$lambda_kappa
      }
    }
    
    Nlambda_beta <- length(lambda_beta)
    Nlambda_kappa <- length(lambda_kappa)
    
    
    # Expand lambda grid:
    lambdas <- expand.grid(kappa = lambda_kappa, beta = lambda_beta)
    Estimates <- vector("list", nrow(lambdas))
    
    ### Algorithm 2 of Rothmana, Levinaa & Ji Zhua
    if (verbose){
      pb <- txtProgressBar(0, nrow(lambdas), style = 3) 
    }
    for (i in seq_len(nrow(lambdas))){
      if ( lambdas$beta[i] == 0 & lambdas$kappa[i] == 0){
        # Unregularized!
        #       SigmaHat <- cov(data, use = "pairwise.complete.obs")
        #       L1Hat <- cov(data_l, data_c, use = "pairwise.complete.obs")
        #       beta <- t(solve(SigmaHat) %*% L1Hat)
        #       kappa <- solve(SigmaHat - beta %*% SigmaHat %*% t(beta))
        X <- data_l
        Y <- data_c
        
        nY <- ncol(Y)
        nX <- ncol(X)
        n <- nrow(X)
        
        beta <- t(Y) %*% X %*% solve(t(X) %*% X)
        
        #####
        ## Compute unconstrained kappa (codes from SparseTSCGM):
        # ZeroIndex <- which(kappa==0, arr.ind=TRUE) ## Select the path of zeros
        S <- 1/(nrow(Y)) * (
          t(Y) %*% Y -
            t(Y) %*% X %*% t(beta) -
            beta %*% t(X) %*% Y +
            beta %*% t(X) %*% X %*% t(beta)
        )
        
        #         out4 <- suppressWarnings(glasso(WS, rho = 0, trace = FALSE))
        # 
        #       kappa <- out4$wi
        S <- (S + t(S)) / 2
        if (any(eigen(S)$value < 0)) stop("Residual covariances not postive definite")
        kappa <- solve(S)
        kappa <- (kappa + t(kappa)) / 2
        
        lik1  = determinant( kappa)$modulus[1]
        lik2 <- sum(diag( kappa%*%S))
        
        pdO = sum(sum(kappa[upper.tri(kappa,diag=FALSE)] !=0))
        pdB = sum(sum(beta !=0))
        
        LLk <-  (n/2)*(lik1-lik2) 
        LLk0 <-  (n/2)*(-lik2)
        
        EBIC <-  -2*LLk + (log(n))*(pdO +pdB) + (pdO  + pdB)*4*gamma*log(2*nY)
        
        #####
        
        Estimates[[i]] <- list(beta = beta, kappa = kappa, EBIC = EBIC)
      } else {
        tryres <- try(Rothmana(data_l, data_c, lambdas$beta[i],lambdas$kappa[i], gamma=gamma,maxit.in=maxit.in, maxit.out = maxit.out,
                               penalize.diagonal = penalize.diagonal,
                               mimic = mimic, likelihood = likelihood)  )
        if (is(tryres,"try-error")){
          Estimates[[i]] <- list(beta=matrix(NA,Nvar,Nvar+1), kappa=matrix(NA,Nvar,Nvar), EBIC = Inf,
                                 error = tryres)
        } else {
          Estimates[[i]] <- tryres
        }
        
      }
      
      if (verbose){
        setTxtProgressBar(pb, i)
      } 
    }
    if (verbose){
      close(pb)
    }
    
    #   
    #   logandbic <- LogLik_and_BIC(data_l, data_c, Estimates)
    #   lambdas$bic <- logandbic$BIC
    #   lambdas$loglik <- logandbic$logLik
    lambdas$ebic <- sapply(Estimates,'[[','EBIC')
    
    if (all(lambdas$ebic==Inf)){
      stop("No model estimated without error")
    }
    # Which minimal BIC:
    min <- which.min(lambdas$ebic)
    Results <- Estimates[[min]]
    
    # Warnings:
    if (length(lambda_beta)>1){
      if (lambdas$beta[[min]] == min(lambda_beta)){
        message("Minimal tuning parameter for beta selected.")
      }
    }
    
    if (length(lambda_kappa)>1){
      if (lambdas$kappa[[min]] == min(lambda_kappa)){
        message("Minimal tuning parameter for kappa selected.")
      }
    }
    
    
    
    # Names of beta:
    colnames(Results$beta) <- colnames(data_l)
    rownames(Results$beta) <- colnames(data_c)
    # Standardize matrices (Wild et al. 2010)
    # partial contemporaneous correlation (PCC) 
    
    Results$PCC <- computePCC(Results$kappa)
    
    # PDS only for lag one!
    if (1 %in% lags){
      Results$PDC <- computePDC(Results$beta[,c("1",paste0(data$vars,"_lag1"))], Results$kappa)  
      if (length(lags) > 1){ 
        warning("Partial directed correlations only computed for lag 1 network.")
      }
    }
    
    Results$path <- lambdas
    Results$labels <- colnames(data_c)
    if (is.null(Results$labels)){
      Results$labels <- paste0("V",seq_len(ncol(data_c)))
    }
    # colnames(Results$beta) <- c("1",Results$labels)
    rownames(Results$beta) <- colnames(Results$kappa) <- rownames(Results$kappa) <-
      colnames(Results$PCC) <- rownames(Results$PCC) <- colnames(Results$PDC) <- rownames(Results$PDC) <-
      Results$labels
    Results$gamma <- gamma
    Results$allResults <- Estimates
    
    Results$N <- nrow(data_c)
    Results$data <- data
    
    class(Results) <- "graphicalVAR"
    
    return(Results)
  }
