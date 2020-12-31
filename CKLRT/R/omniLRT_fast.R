#' Composite kernel machine regression based likelihood ratio test.
#'
#' Composite kernel machine regression based likelihood ratio test. The approximate method for likelihood ratio test tend to be too conservative for small alpha values. We recommend not using it in GWAS
#' @param y : y is the vecgtor of the continous outcomes.
#' @param X : X denotes the additional covariates.
#' @param K1 : K1 is the first kernel corresponding to the genetic main effect.
#' @param K2 : K2 is the second kernel corresponding to the genetic and environment interaction effect.
#' @param N : N is th total number of randomly generated normal variables used to generate the emprical null distribution of LRT. Default value is 10,000.
#' @param length.lambda : the length of lambda. Dafult value is 200. The values of lambda are all more than 0.
#' @param length.rho : the length of rho. Default value is 21. The values of rho are between 0 and 1.
#'
#' @return the result is a list containing three elements. 1. p.dir is the p-value of likelihood ratio test based on emprical distrition. 2. p.aud is the p-value by approximating the null distribution as a mixture of a point mass at zero with probability b and weighted chi square distribution with d degrees of freedom with probality of 1-b. 3. LR is the likelihood ratio test statistics.
#' @export
#'
#' @examples
#' set.seed(6)
#' n = 50 # the number of observations
#' X = rnorm(n) # the other covariates
#' p = 2 # two snp in a gene will be simulated
#' G = runif(n*p)< 0.5
#' G = G + runif(n*p) < 0.5
#' G = matrix(G, n,p) #genetic matrix
#' E = (runif(n) < 0.5)^2 #enviroment effect
#' y = rnorm(n) + G[,1] * 0.3 #observations
#' omniLRT_fast(y, X =  cbind(X, E),K1 = G %*% t(G),K2 = (G*E) %*% t(G * E))
#' @importFrom MASS ginv
#' @import nlme
#' @import Rcpp
#' @import mgcv
#' @import compiler
#' @import stats



omniLRT_fast = function(y, X,K1, K2, N = 10000,length.lambda = 200, length.rho = 21){
  method = "ML"
  Lambdas = exp(seq(from = -12, to = 12, length.out = length.lambda))
  all_rho = seq(from = 0,to = 1, length.out = length.rho)
  n  = length(y)
  if (is.null(X)){
    X1 = matrix(1, nrow=n)
    px = 1
  }else{
    X1 = cbind(1,X)
    px = ncol(X1)
  }
  XX   = MatMult_C(t(X1),X1)
  
  P0   = diag(n)- MatMult_C(MatMult_C(X1,ginv(XX)),t(X1))
  
  eP   = Eigen_C(P0)
  A    = eP$vector[,eP$values > 1e-10]
  # invP = ginv(P0)
  # A2   = eP$vectors %*% diag(eP$values)
  eK1  = Eigen_C(K1)
  wK1  = which(eK1$values  > 1e-10)
  # phi1 = t(t(eK1$vectors[,wK1])*sqrt(eK1$values[wK1]))
  if (length(wK1) == 1){
    phi1 = eK1$vectors[,wK1] * sqrt(eK1$values[wK1])
  }else{
    phi1 = t(t(eK1$vectors[,wK1])*sqrt(eK1$values[wK1]))
  }
  
  eK2  = Eigen_C(K2)
  wK2  = which(eK2$values  > 1e-10)
  
  if (length(wK2) == 1){
    phi2 = eK2$vectors[,wK2] * sqrt(eK2$values[wK2])
  }else{
    
    phi2 = t(t(eK2$vectors[,wK2])*sqrt(eK2$values[wK2]))
  }
  
  group= rep(1,n)
  if(is.null(X)){
    fit1 = lme(y~1, random = list(group=pdIdent(~-1+phi1), group = pdIdent(~-1+phi2)), method = "ML") #  Default = REML
    fit0 = stats::lm(y~1)
    LR = max(0, 2*(stats::logLik(fit1, REML = F) -stats::logLik(fit0, REML = F)))
    
  }else{
    fit1 = lme(y~X, random = list(group=pdIdent(~-1+phi1), group = pdIdent(~-1+phi2)), method = "ML") #  Default = REML
    fit0 = stats::lm(y~X)
    LR = max(0, 2*(stats::logLik(fit1, REML = F) -stats::logLik(fit0, REML = F)))
    
  }
  
  if (LR <= 0){
    p.dir = 1;
    p.au1=1; p.aud = 1
    
  }else{
    #For the first kernel
    LR0_allRho = matrix(NA, N, length.rho)
    #set.seed(123)
    w = matrix(stats::rnorm(N*(n-px)), n-px,N)
    LR0_fixRho = matrix(NA, N, length.lambda)
    # The baseline W1
    rho = 0   # K = rho*K1 + (1-rho)*K2 = K2
    K = K2
    k = length(wK2)
    xi = eK2$values[wK2]
    
    mu = Eigen_C_value(t(phi2) %*% P0%*% phi2)
    # mu = eigen(P0 %*% K2)$values
    
    mu = mu/max(mu, xi)
    xi = xi/max(mu,xi)
    AKA = MatMult_C(MatMult_C(t(A),K),A)
    eV = Eigen_C(AKA)
    U_1 = eV$vectors
    w.double = w^2
    w1 = (w.double)[1:k,]
    w2 = ColSum_C((w.double)[-(1:k),])
    
    if (length(mu) < k){mu = c(mu,rep(0, k - length(mu)))}
    if (length(xi) < k){xi = c(xi,rep(0, k - length(xi)))}
    
    LR0_fixRho <- LR0_fixRho_LRT_C(Lambdas,
                                   mu,
                                   w1,
                                   w2,
                                   n-px,
                                   xi)
    LR0_allRho[,1] = MatrixRowMax_C(LR0_fixRho)
    LR0_allRho <- doubleloop_LRT(K1,
                                 K2,
                                 P0,
                                 A,
                                 U_1,
                                 w,
                                 Lambdas,
                                 n-px,
                                 all_rho,
                                 LR0_allRho)
    LR0 = MatrixRowMax_C(LR0_allRho)
    LR0 = ifelse(LR0 > 0, LR0, 0)
    p.dir = mean(LR < LR0)
    p.aud= getp_aud_estimate_pi_first(null = LR0, LR = LR)$p
    
  }
  out = list(p.dir = p.dir,p.aud = p.aud, LR = LR)
  return(out)
}


