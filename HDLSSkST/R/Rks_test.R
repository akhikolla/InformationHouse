# -----------------------------------rand_index--------------------------------------
randfun <- function(lvel, dv) {
  lem <- abs(sapply(lvel, function(x) x-lvel))
  lem[lem > 1] <- 1
  olem <- abs(sapply(dv, function(x) x-dv))
  olem[olem > 1] <- 1
  ada <- sum(abs(lem-olem))/2
  tnp <- choose(dim(lem)[1],2)
  ada/tnp
}
# ------------------------------------------------------------------------------------
#----------------------------pmf of G-hyg distribution--------------------------------
pmf <- function(M){
  cols <- colSums(M)
  nume <- prod(as.numeric(unlist(apply(M, 1, function(v) dmultinom(v, size = NULL, rep(1/ncol(M),ncol(M)), log = FALSE)))))
  deno <- dmultinom(cols, size = NULL, rep(1/ncol(M),ncol(M)), log = FALSE)
  nume/deno  
}
#-------------------------------------------------------------------------------------


# ------------------------------------rand_index_test-----------------------------------------------------------
RItest <- function(M, labels, sizes, randomization = TRUE, clust_alg = "knwClustNo", s_fn = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  if(is.data.frame(M))
    M <- as.matrix(M)
  if(!is.numeric(M) || !is.numeric(labels) || !is.numeric(sizes) || anyNA(M) || anyNA(labels) || anyNA(sizes))
    stop("all entries of obsevations, cluster membership index of observations and number of observations from each of the population must be finite")
  if(any(sizes < 2L))
    stop("Each of population must have at least 2 obsevations")
  if(nrow(M)!=length(labels)||length(labels)!=sum(sizes))
    stop("Number of observations and number of membership of observations must be same")
  N <- nrow(M)
  # d <- ncol(M)
  n_clust <- length(sizes)
  if(clust_alg=="knwClustNo"){
    dvec <- as.numeric(gMADD(s_fn, n_clust, lb, M))
  }
  else if(clust_alg=="estClustNo"){
    maxClNo <- 2*n_clust
    dvec_di_mat <- as.matrix(gMADD_DI(s_fn, maxClNo, lb, M)) 
    est_no_cl <- which.max(dvec_di_mat[ ,(N+1)])
    dvec <- dvec_di_mat[est_no_cl,1:N]
  }
  obs_ct <- unname(table(labels, dvec))
  rand_id <- randfun(labels, dvec)

  sim_ri <- rep(0, n_sts)
  for(j in 1:n_sts){
    sc_tab <- as.matrix(rctab(obs_ct))
    s_dvec <- as.numeric(c(unlist(apply(sc_tab, 1, function(v) rep(0:(ncol(sc_tab)-1), v)))))
    sim_ri[j] <- randfun(labels, s_dvec)
  }
  
  pvalue_RI <-  mean(sim_ri<=rand_id)
  # R_alpha <- as.numeric(quantile(sim_ri, alpha))
  R_alpha <- sort(sim_ri)[floor(alpha*n_sts)]
  dec_r_x <- 0
  # if(rand_id <= R_alpha ){dec_r_x <- 1}
  
  ri_gamma <- 0
 if(randomization==TRUE){
   
   if(rand_id < R_alpha){
     dec_r_x <- 1
   } 
   
   if(rand_id==R_alpha){
     u <-runif(1) 
     s0 <- mean(sim_ri==R_alpha)
     s1 <- mean(sim_ri<R_alpha)
     
     ri_gamma <- (alpha - s1)/s0     
     if(u < ri_gamma){dec_r_x <- 1}
   }
 }
 else if(randomization==FALSE){
   if(rand_id < R_alpha){
     dec_r_x <- 1
   } 
 }
  #output
  output<-list()
  if(clust_alg=="knwClustNo"){
    output$estClustLabel <- dvec
    output$obsCtyTab <- obs_ct
    output$ObservedRI <- rand_id
    output$RICutoff <- R_alpha
    output$randomGamma <- ri_gamma
    output$estPvalue <- pvalue_RI
    output$decisionRI <- dec_r_x
  }
  else if(clust_alg=="estClustNo"){
    output$estClustNo <- est_no_cl
    output$estClustLabel <- dvec
    output$obsCtyTab <- obs_ct
    output$ObservedRI <- rand_id
    output$RICutoff <- R_alpha
    output$randomGamma <- ri_gamma
    output$estPvalue <- pvalue_RI
    output$decisionRI <- dec_r_x
  }
  return(output)
}
# ------------------------------------------------------------------------------------------------

# ------------------------------Fisher_Exact_Independence_test------------------------------------------------------
FEItest <- function(M, labels, sizes, randomization = TRUE, clust_alg = "knwClustNo", s_fn = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  if(is.data.frame(M))
    M <- as.matrix(M)
  if(!is.numeric(M) || !is.numeric(labels) || !is.numeric(sizes) || anyNA(M) || anyNA(labels) || anyNA(sizes))
    stop("all entries of obsevations, cluster membership index of observations and number of observations from each of the population must be finite")
  if(any(sizes < 2L))
    stop("Each of population must have at least 2 obsevations")
  if(nrow(M)!=length(labels)||length(labels)!=sum(sizes))
    stop("Number of observations and number of membership of observations must be same")
  N <- nrow(M)
  # d <- ncol(M)
  n_clust <- length(sizes)
  if(clust_alg=="knwClustNo"){
    dvec <- as.numeric(gMADD(s_fn, n_clust, lb, M))
  }
  else if(clust_alg=="estClustNo"){
    maxClNo <- 2*n_clust
    dvec_di_mat <- as.matrix(gMADD_DI(s_fn, maxClNo, lb, M)) 
    est_no_cl <- which.max(dvec_di_mat[ ,(N+1)])
    dvec <- dvec_di_mat[est_no_cl,1:N]
  }
  obs_ct <- unname(table(labels, dvec))
  fpmf <- pmf(obs_ct)

  fsim_pmf <- rep(0, n_sts)
  for(j in 1:n_sts){
    fsc_tab <- as.matrix(rctab(obs_ct))
    fsim_pmf[j] <- pmf(fsc_tab)
  }
  
  pvalue_FEI <- mean(fsim_pmf<=fpmf)
  # f_alpha <- as.numeric(quantile(fsim_pmf, alpha))
  f_alpha <- sort(fsim_pmf)[floor(alpha*n_sts)]
  decfisher <- 0
  # if(fpmf <= f_alpha){decfisher <- 1}
  
  fgamma <- 0
  if(randomization==TRUE){
    
    if(fpmf < f_alpha){
      decfisher <- 1
    } 
    
    if(fpmf==f_alpha){
      u1 <-runif(1) 
      fs0 <- mean(fsim_pmf==f_alpha)
      fs1 <- mean(fsim_pmf<f_alpha)
      
      fgamma <- (alpha - fs1)/fs0     
      if(u1 < fgamma){decfisher <- 1}
    }
  }
  else if(randomization==FALSE){
    if(fpmf < f_alpha){
      decfisher <- 1
    } 
  }
  #output
  output<-list()
  if(clust_alg=="knwClustNo"){
    output$estClustLabel <- dvec
    output$obsCtyTab <- obs_ct
    output$ObservedProb <- fpmf
    output$FCutoff <- f_alpha
    output$randomGamma <- fgamma
    output$estPvalue <- pvalue_FEI
    output$decisionF <- decfisher
  }
  else if(clust_alg=="estClustNo"){
    output$estClustNo <- est_no_cl
    output$estClustLabel <- dvec
    output$obsCtyTab <- obs_ct
    output$ObservedProb <- fpmf
    output$FCutoff <- f_alpha
    output$randomGamma <- fgamma
    output$estPvalue <- pvalue_FEI
    output$decisionF <- decfisher
  }
  return(output)
}

# ------------------------------------------------------------------------------------------

# ------------------------Aggregate_rand_index_test----------------------------------------
ARItest <- function(M, sizes, randomization = TRUE, clust_alg = "knwClustNo", s_fn = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  if(is.data.frame(M))
    M <- as.matrix(M)
  if(!is.numeric(M) || !is.numeric(sizes) || anyNA(M) || anyNA(sizes))
    stop("all entries of obsevations and number of observations from each of the population must be finite")
  if(any(sizes < 2L))
    stop("Each of population must have at least 2 obsevations")
  no_clust <- length(sizes)
  cumS <- rep(0,no_clust+1)
  cumS[1] <- 0
  for(kk in 2:(no_clust+1)){
    cumS[kk] <- cumS[kk-1] + sizes[kk-1]
  }
  
  aggRI <- 1
  epvalues <- 0
  minrow <- rep(1,n_sts)
  for(i in 1:(no_clust-1)){
    mm <- sizes[i]
    for(j in (i+1):no_clust){
      nn <- sizes[j]
      level <- c(rep(0,mm), rep(1,nn))
      obsM <- M[c((cumS[i]+1):cumS[i+1],(cumS[j]+1):cumS[j+1]), ]
      testFun <- RItest(obsM, labels=level, sizes = c(mm,nn), randomization, clust_alg, s_fn, lb, n_sts, alpha) 
      aggRI <- min(aggRI,testFun$ObservedRI)
      epvalues <- c(epvalues,testFun$estPvalue)
      
      sim_ri <- rep(0, n_sts)
      for(jj in 1:n_sts){
        sc_tab <- as.matrix(rctab(testFun$obsCtyTab))
        s_dvec <- c(unlist(apply(sc_tab, 1, function(v) rep(0:(ncol(sc_tab)-1), v))))
        sim_ri[jj] = randfun(level, s_dvec)
      }
      minrow <- pmin(minrow,sim_ri)
    }
  }
  
  # AR_alpha <- as.numeric(quantile(minrow, alpha))
  AR_alpha <- sort(minrow)[floor(alpha*n_sts)]
  dec_ARI <- 0
  # if(aggRI <= AR_alpha ){dec_ARI <- 1}
  
  randGamma <- 0
  if(randomization==TRUE){
    
    if(aggRI< AR_alpha){
      dec_ARI <- 1
    } 
    
    if(aggRI==AR_alpha){
      u <-runif(1)
      s0 <- mean(minrow==AR_alpha)
      s1 <- mean(minrow<AR_alpha)
      
      randGamma <- (alpha - s1)/s0     
      if(u < randGamma){dec_ARI <- 1}
    }
  }
  else if(randomization==FALSE){
    if(aggRI < AR_alpha){
      dec_ARI <- 1
    } 
  }
  lpvs <- no_clust*(no_clust-1)/2
  pvalues <- epvalues[2:(lpvs+1)]
  criticalValues <- sapply(1:lpvs, function(i) alpha/(lpvs-i+1))
  notequal <-  (sort(pvalues) < criticalValues)
  orp <- as.vector(order(pvalues))
  idx <- as.vector(which(notequal))
  notequalpop <- orp[idx]
  hyposet <- t(combn(1:no_clust,2))
  multipletest <- data.frame(Population = hyposet[notequalpop, ], pvalues=pvalues[notequalpop], criticalValues=criticalValues[notequalpop])
  
  #output
  output<-list()
  output$ARIStat <- aggRI
  output$ARICutoff <- AR_alpha
  output$randomGamma <- randGamma
  output$decisionARI <- dec_ARI
  output$multipleTest <- multipletest
 
  # if(clust_alg=="knwClustNo"){
  #   output$ARIStat <- aggRI
  #   output$ARICutoff <- AR_alpha
  #   output$randomGamma <- randGamma
  #   output$decisionARI <- dec_ARI
  #   output$multipleTest <- multipletest
  # }
  # else if(clust_alg=="estClustNo"){
  #   output$ARIStat <- aggRI
  #   output$ARICutoff <- AR_alpha
  #   output$randomGamma <- randGamma
  #   output$decisionARI <- dec_ARI
  #   output$multipleTest <- multipletest
  #   }
  return(output)
}
  
# ------------------------------------------------------------------------------------------

# ------------------------Aggregate_Fisher_Exact_Independence_test------------------------------
AFEItest <- function(M, sizes, randomization = TRUE, clust_alg = "knwClustNo", s_fn = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  if(is.data.frame(M))
    M <- as.matrix(M)
  if(!is.numeric(M) || !is.numeric(sizes) || anyNA(M) || anyNA(sizes))
    stop("all entries of obsevations and number of observations from each of the population must be finite")
  if(any(sizes < 2L))
    stop("Each of population must have at least 2 obsevations")
  no_clust <- length(sizes)
  cumS <- rep(0,no_clust+1)
  cumS[1] <- 0
  for(kk in 2:(no_clust+1)){
    cumS[kk] <- cumS[kk-1] + sizes[kk-1]
  }
  
  aggAFEI <- 1
  epvalues <- 0
  minrow <- rep(1,n_sts)
  
  for(i in 1:(no_clust-1)){
    mm <- sizes[i]
    for(j in (i+1):no_clust){
      nn <- sizes[j]
      level <- c(rep(0,mm), rep(1,nn))
      obsM <- M[c((cumS[i]+1):cumS[i+1],(cumS[j]+1):cumS[j+1]), ]
      testFun <- FEItest(obsM, labels=level, sizes = c(mm,nn), randomization, clust_alg, s_fn, lb, n_sts, alpha) 
      aggAFEI <- min(aggAFEI,testFun$ObservedProb)
      epvalues <- c(epvalues,testFun$estPvalue)
      
      fsim_pmf <- rep(0, n_sts)
      for(jj in 1:n_sts){
        fsc_tab <- as.matrix(rctab(testFun$obsCtyTab))
        fsim_pmf[jj] = pmf(fsc_tab)
      }
      
      minrow <- pmin(minrow,fsim_pmf)
    }
  }
  
  # AF_alpha <- as.numeric(quantile(minrow, alpha))
  AF_alpha <- sort(minrow)[floor(alpha*n_sts)]
  dec_AFEI <- 0
  # if(aggAFEI <= AF_alpha ){dec_AFEI <- 1}
  
  randGamma <- 0
  if(randomization==TRUE){
    
    if(aggAFEI< AF_alpha){
      dec_AFEI <- 1
    } 
    
    if(aggAFEI==AF_alpha){
      u1 <-runif(1) 
      fs0 <- mean(minrow==AF_alpha)
      fs1 <- mean(minrow<AF_alpha)
      
      randGamma <- (alpha - fs1)/fs0     
      if(u1 < randGamma){dec_AFEI <- 1}
    }
  }
  else if(randomization==FALSE){
    if(aggAFEI< AF_alpha){
      dec_AFEI <- 1
    } 
  }
  lpvs <- no_clust*(no_clust-1)/2
  pvalues <- epvalues[2:(lpvs+1)]
  criticalValues <- sapply(1:lpvs, function(i) alpha/(lpvs-i+1))
  notequal <-  (sort(pvalues) < criticalValues)
  orp <- as.vector(order(pvalues))
  idx <- as.vector(which(notequal))
  notequalpop <- orp[idx]
  hyposet <- t(combn(1:no_clust,2))
  multipletest <- data.frame(Population = hyposet[notequalpop, ], pvalues=pvalues[notequalpop], criticalValues=criticalValues[notequalpop])
  
  #output
  output<-list()
  output$AFEIStat <- aggAFEI
  output$AFCutoff <- AF_alpha
  output$randomGamma <- randGamma
  output$decisionAFEI <- dec_AFEI
  output$multipleTest <- multipletest
  
  # if(clust_alg=="knwClustNo"){
  #   output$AFEIStat <- aggAFEI
  #   output$AFCutoff <- AF_alpha
  #   output$randomGamma<- randGamma
  #   output$decisionAFEI <- dec_AFEI
  #   output$multipleTest <- multipletest
  # }
  # else if(clust_alg=="estClustNo"){
  #   output$AFEIStat <- aggAFEI
  #   output$AFCutoff <- AF_alpha
  #   output$randomGamma <- randGamma
  #   output$decisionAFEI <- dec_AFEI
  #   output$multipleTest <- multipletest
  #   }
  return(output)
}
# ------------------------------------------------------------------------------------------

