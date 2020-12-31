#Utility functions for the DriftBurstHypothesis R package

### polyfit2d
fit2d3rdDegree = function(A, X, Y){
  
  X = matrix(rep(as.numeric(X),4), ncol = 4)
  Y = matrix(Y, ncol = 16, nrow = length(Y), byrow = FALSE)
  Z = matrix(1, ncol = ncol(Y), nrow = nrow(Y))
  
  for (i in 1:3) {
    Z[,1:i] = Z[,1:i] * X[,1:i]
  }
  
  for (i in 1:3) {
    
    foo = 4*i
    Z[, 1:foo] = Z[, 1:foo] * Y[, 1:foo]
    for (j in 1:3) {
      Z[, (foo+1):(foo+j)] = Z[, (foo+1):(foo+j)] * X[,1:j]
    }
  }
  p = mldivide(Z,as.numeric(A))
  return(p)
}






DBHCriticalValues = function(x, alpha){
  tStat = getDB(x)
  alpha_used = alpha
  
  nObs = length(tStat)
  rho = DBHSysData$rho
  alpha = DBHSysData$alpha
  m = DBHSysData$m
  arr = array(DBHSysData$arr, dim = c(23, 14, 6), dimnames = list('m' = m, 'rho' = rho,'alpha' = alpha))  
  
  rho_used = acf(tStat, plot = FALSE)$acf[2]
  A = matrix(NA, ncol = nrow(arr), nrow = nrow(arr))
  P = matrix(NA, ncol = dim(arr)[3], nrow = 16)
  rho = as.numeric(dimnames(arr)$rho)
  add_rho = c(0.75,0.85,0.91,0.92,0.93,0.94,0.96,0.97,0.98)
  rho_i = sort(c(rho, add_rho))
  logRho = log(1-rho)
  logRho_i = log(1-rho_i)
  logM = matrix(log(as.numeric(dimnames(arr)$m)), ncol = nrow(arr), nrow = nrow(arr), byrow = TRUE)
  logRho_mat = matrix(logRho_i, ncol = nrow(arr), nrow = nrow(arr), byrow = FALSE)
  
  for (i in 1:length(alpha)) {
    
    for (j in 1:nrow(arr)) {
      A[ ,j] = approx(logRho, arr[j, ,i], logRho_i)$y
    }
    P[,i] = fit2d3rdDegree(A, logM, logRho_mat)
  }
  
  idx = alpha %in% alpha_used
  normalizedQuantile = P[1,idx]
  for (i in 1:3) {
    normalizedQuantile = normalizedQuantile * log(nObs) + P[1+i, idx]
  }
  
  
  for (i in 1:3) {
    foo = 4 * i+1
    g = P[foo, idx]
    
    for (j in 1:3) {
      g = g * log(nObs) + P[j+ foo, idx]
    }
    normalizedQuantile = normalizedQuantile * log(1 - rho_used) + g
    
  }
  
  AM = sqrt(2 * log(nObs))
  BM = AM - 0.5 * log(pi * log(nObs)) / AM
  
  quantile = normalizedQuantile / AM + BM
  
  return(list("normalizedQuantile" = normalizedQuantile, "quantile" = quantile))
}



