#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#define _USE_MATH_DEFINES
#include <math.h>


NumericMatrix duplicateMat(NumericMatrix X){
  NumericMatrix Y(X.nrow(),X.ncol());
  for (int i=0; i<X.nrow(); i++){
    for (int j=0; j < X.ncol(); j++){
      Y(i,j) = X(i,j);
    }
  }
  return(Y);
}

double signfun(double x){
  double res = 0;
  if (x < 0){
    res = -1;
  }
  if (x > 0){
    res = 1;
  }  
  return(res);
}

double noneg(double x){
  if (x < 0){
    x = 0;
  }
  return(x);
}


double absfun(double x){
  if (x < 0){
    x = -1*x;
  }
  return(x);
}

// Ridge regression beta
// [[Rcpp::export]]
NumericMatrix beta_ridge_C(NumericMatrix X, NumericMatrix Y, double lambda_beta){
  
  int n = X.nrow(), nX = X.ncol(), nY = Y.ncol();
  // Create lambda identity:
  NumericMatrix lambdaIden(nX,nX);
  std::fill(lambdaIden.begin(), lambdaIden.end(), 0.0);
  for (int i=0;i<nX;i++){
    lambdaIden(i,i) = lambda_beta;
  }
  
  arma::mat X_(X.begin(), n, nX, false);       // reuses memory and avoids extra copy
  arma::mat Y_(Y.begin(), n, nY, false);       // reuses memory and avoids extra copy
  arma::mat lambdaIden_(lambdaIden.begin(), nX, nX, false);       // reuses memory and avoids extra copy
  
  // Compute ridge beta:
  arma::mat beta_ridge_ = inv(trans(X_) * X_ + lambdaIden_) * trans(X_) * Y_;
  
  return(Rcpp::as<Rcpp::NumericMatrix>(wrap(beta_ridge_)));
}


// Compute Beta given Kappa:
// [[Rcpp::export]]
NumericMatrix Beta_C(NumericMatrix kappa, NumericMatrix beta, NumericMatrix X, NumericMatrix Y, 
double lambda_beta, NumericMatrix lambda_beta_mat, double convergence, int maxit){
  
  int n = X.nrow(), nX = X.ncol(), nY = Y.ncol();
  
  // Convert matrices without reusing memory:
  arma::mat X_(X.begin(), n, nX, false);       // reuses memory and avoids extra copy
  arma::mat Y_(Y.begin(), n, nY, false);       // reuses memory and avoids extra copy
  arma::mat kappa_(kappa.begin(), nY, nY, false);       // reuses memory and avoids extra copy
  
  // new matrices:
  arma::mat S_ = trans(X_) * X_;
  NumericMatrix S = as<Rcpp::NumericMatrix>(wrap(S_));
  arma::mat H_ = trans(X_) * Y_ * kappa_;
  NumericMatrix H = as<Rcpp::NumericMatrix>(wrap(H_));
  
  
  // Ridge:
  NumericMatrix beta_ridge = beta_ridge_C(X, Y, lambda_beta);
  double ridgecriterium = 0;
  for (int j=0; j<beta.nrow(); j++){
    for (int k=0; k<beta.ncol(); k++){
      ridgecriterium += absfun(beta_ridge(j,k));
    }
  }
  double criterium;
  
  NumericMatrix beta_new = duplicateMat(beta);
  int it = 0;
  // Store beta:
  do{
    NumericMatrix beta_old = duplicateMat(beta_new);
    
    // Sequential update:
    for (int r=0; r<beta_new.nrow();r++){
      for (int c=0; c<beta_new.ncol();c++){
        double u = 0;
        for (int j=0;j<beta_new.nrow();j++){
          for (int k=0; k<beta_new.ncol(); k++){
            u += beta_new(j, k) * S(r, j) * kappa(k, c);
          }
        }
        beta_new(r,c) = signfun(beta_new(r,c) + (H(r,c) - u)/(S(r,r) * kappa(c,c))) * noneg(absfun(beta_new(r,c) + (H(r,c) - u)/(S(r,r)*kappa(c,c))) - n*lambda_beta_mat(r,c)/(S(r,r)*kappa(c,c)));
      }
    }
    
    criterium = 0;
    for (int j=0; j<beta_new.nrow(); j++){
      for (int k=0; k<beta_new.ncol(); k++){
        criterium += absfun(beta_new(j,k) - beta_old(j,k));
      }
    }
    
    it++;
  } while (it < maxit && criterium > (convergence * ridgecriterium));
  
  if (it >= maxit){
    Rcpp::Rcout << "\nModel did NOT converge in inner loop";
  }
  
  return(beta_new);
}



// [[Rcpp::export]]
double VAR_logLik_C(NumericMatrix X, NumericMatrix Y, NumericMatrix kappa, NumericMatrix beta){
  // http://webspace.qmul.ac.uk/aferreira/lect2-var2_handout.pdf
  
  int T = X.nrow();
  int nX = X.ncol(), nY = Y.ncol();
  
  // Convert matrices without reusing memory:
  arma::mat X_arma(X.begin(), T, nX, false);      
  arma::mat Y_arma(Y.begin(), T, nY, false);      
  arma::mat kappa_arma(kappa.begin(), nY, nY, false);     
  arma::mat beta_arma(beta.begin(), nY, nX, false);   
  
  
  
  double sum = 0;
  for (int t=0; t<T; t++){
    arma::mat foo = trans(trans(Y_arma.row(t)) - beta_arma * trans(X_arma.row(t))) * kappa_arma * (trans(Y_arma.row(t)) - beta_arma * trans(X_arma.row(t)));
    sum += foo[0];
  }
  
  double res =  -((double)T * (double)nY/2) * log(2*M_PI) + ((double)T/2) * log(det(kappa_arma)) - 0.5 * sum;
  
  return(res);
}

// [[Rcpp::export]]
List LogLik_and_BIC(NumericMatrix X, NumericMatrix Y, List estimates){

  int n = X.nrow();
  
  int N = estimates.size();
  NumericVector LogLiks(N);
  NumericVector BICs(N);
  
  for (int k=0; k<N; k++){
    
    List el = estimates[k];
    NumericMatrix kappa = el["kappa"];
    NumericMatrix beta = el["beta"];

    LogLiks[k] = VAR_logLik_C(X, Y, kappa, beta);
    
    // Number of parameters:
    int nPar = 0;
    for (int i=0; i<kappa.nrow(); i++){
      for (int j=i; j<kappa.ncol(); j++){
        if (i != j){
          if (kappa(i,j) != 0){
            nPar++;
          }
        }
      }
    }
    for (int i=0; i<beta.nrow(); i++){
      for (int j=0; j<beta.ncol(); j++){
        if (beta(i,j) != 0){
          nPar++;
        }
      }
    }
    
    BICs[k] = -2 * LogLiks[k] + nPar * log((double)n);
  }
  
  List Results;
  Results["logLik"] = LogLiks;
  Results["BIC"] = BICs;
  
  return(Results);
}



/*

# Ridge estimate:
#   beta_ridge <- matrix(NA, Nvar, Nvar)
#   for (i in seq_len(Nvar)){
#     glmres <- glmnet(X, Y[,i], alpha = 0, lambda = lambda_beta)
#     beta_ridge[i,] <- coef(glmres)[-1]
#   }
beta_ridge <- solve(t(X)%*%X + lambda_beta*diag(Nvar))%*%t(X)%*%Y

repeat{
beta_old <- beta

# Random sequence of row and column:
seq <- expand.grid(r=seq_len(Nvar), c=seq_len(Nvar))
seq <- seq[sample(seq_len(nrow(seq))),]

for (i in seq_len(nrow(seq))){
r <- seq$r[i]
c <- seq$c[i]

u <- sum(t(beta) * kappa[,c,drop=FALSE] %*% S[r,,drop=FALSE] )

beta[r,c] <- sign(beta[r,c] + (H[r,c] - u)/(S[r,r]*kappa[c,c])) * max(abs(beta[r,c] + (H[r,c] - u)/(S[r,r]*kappa[c,c])) - n*lambda_beta/(S[r,r]*kappa[c,c]),0)
}

if (sum(abs(beta - beta_old)) < (convergence * sum(abs(beta_ridge)))){
break
}
}

return(beta)
}


*/
