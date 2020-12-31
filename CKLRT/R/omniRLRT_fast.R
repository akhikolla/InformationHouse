#' Composite kernel machine regression based restricted likelihood ratio test
#'
#' @param y vector of the continous outcomes.
#' @param X the additional covariates.
#' @param K1 the first kernel corresponding to the genetic main effect.
#' @param K2 the second kernel corresponding to the genetic and environment interaction effect.
#' @param N total number of randomly generated normal variables used to generate the emprical null distribution of LRT. Default value is 10,000.
#' @param length.lambda the length of lambda. Dafult value is 200. The values of lambda are all more than 0.
#' @param length.rho the length of rho. Default value is 21. The values of rho are between 0 and 1.
#'
#' @return the result is a list containing three elements. 1. p.dir is the p-value of restricted likelihood ratio test based on emprical distrition. 2. p.aud is the p-value by approximating the null distribution as a mixture of a point mass at zero with probability b and weighted chi square distribution with d degrees of freedom with probality of 1-b. 3. LR is the likelihood ratio test statistics.

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
#' omniRLRT_fast(y, X =  cbind(X, E),K1 = G %*% t(G),K2 = (G*E) %*% t(G * E))
#' @importFrom MASS ginv
#' @import nlme
#' @import Rcpp
#' @import mgcv
#' @import compiler

omniRLRT_fast = function(y, X,K1, K2, N = 10000, length.rho = 200, length.lambda = 21){
  method = "REML"
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
    phi1 = t(t(eK1$vectors[,wK1])*sqrt(eK1$values[wK1])) # this actually works for all of them
  }

  eK2  = Eigen_C(K2)
  wK2  = which(eK2$values  > 1e-10)

  if (length(wK2) == 1){
    phi2 = eK2$vectors[,wK2] * sqrt(eK2$values[wK2])
  }else{
    phi2 = t(t(eK2$vectors[,wK2])*sqrt(eK2$values[wK2]))
  }
  # if wK1 and wK2 are 1, what would this happen to others?

  group= rep(1,n)
  fit1 = lme(y~X, random = list(group=pdIdent(~-1+phi1), group = pdIdent(~-1+phi2))) #  Default = REML
  fit0 = lm(y~X)
  LR = max(0, 2*(logLik(fit1, REML = T) -logLik(fit0, REML = T)))

  if (LR <= 0){
    p.dir =p.au1= p.aud = 1
  }else{

    #For the first kernel
    LR0_allRho = matrix(NA, N, length.rho)
    #set.seed(123)
    w = matrix(rnorm(N*(n-px)), n-px,N)
    LR0_fixRho = matrix(NA, N, length.lambda)

    rho = 0
    K = K2
    k = length(wK2)
    xi = eK2$values[wK2]
    AKA = MatMult_C(MatMult_C(t(A),K),A)
    eV = Eigen_C(AKA)
    #eV = Eigen_C(AKA)
    U_1 = eV$vectors
    mu = eV$values[eV$values > 1e-10]


    mu = mu/max(mu, xi)
    xi = xi/max(mu,xi)

    # original U_1 in this case
    W1 = MatMult_C(A,U_1)
    w.double <- w^2
    w1 = w.double[1:k,]
    w2 = ColSum_C((w.double)[-(1:k),])

    if (length(mu) < k){mu = c(mu,rep(0, k - length(mu)))}
    if (length(xi) < k){xi = c(xi,rep(0, k - length(xi)))}


    LR0_fixRho <- LR0_fixRho_C(Lambdas,
                               mu,
                               w1,
                               w2,
                               n-px)
    # for (i in 1:length.lambda){
    #   lam = Lambdas[i]
    #   Dn = (1/(1 + lam*mu))%*%w1+ w2
    #   Nn = (lam*mu/(1 + lam*mu))%*%w1
    #   temp = (n-px)*log(1 + Nn/Dn) - Sum_C(log(1 + lam*mu))
    #   LR0_fixRho[,i] = ifelse(temp < 0, 0, temp)
    # }
    LR0_allRho[,1] = MatrixRowMax_C(LR0_fixRho)

    LR0_allRho <- doubleloop(K1,
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
    p.au1 = getp_au1(null = LR0, LR = LR)$p
    p.aud= getp_aud_estimate_pi_first(null = LR0, LR = LR)$p
  }
  out = list(p.dir = p.dir,p.aud = p.aud, LR = LR)
  return(out)
  }


