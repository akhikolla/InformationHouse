/******************  created by Erlis Ruli 14/01/2014  ************************
****************** Robust Approximate Bayesian Infernece  *********************
******************************************************************************/

#include "extrafuns.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

// Indipendent Normal-Half Cauchy
//[[Rcpp::export]]
double
  dPrior_lmm(arma::vec betas,
             double log_sig2_b,
             double log_sig2_eps,
             double beta_stDev,
             double sig2_scale,
             int n_betas,
             bool log_){
  double ans = 0.0;
  for(int i=0; i<n_betas; i++){
    ans += Rf_dnorm4(betas(i), 0.0, beta_stDev, log_);
  }
  ans +=  dhalfCauchy(exp(log_sig2_b), sig2_scale, log_) + log_sig2_b +
          dhalfCauchy(exp(log_sig2_eps), sig2_scale, log_) + log_sig2_eps;

  if(log_ == false) {
    ans = exp(ans);
  }

  return ans;
}

// product of dim_x normal densities with specified variance h
//[[Rcpp::export]]
double
dmvnorm_I(arma::vec x, double h, int dim_x, bool log_){
  double ans = 0.0;
  for(int i = 0; i< dim_x; i++){
    ans += Rf_dnorm4(x(i), 0.0, h, log_);
  }
  if(log_ == false) {
    ans = exp(ans);
  }
  return ans;
}

// V matrix and related quantites
//[[Rcpp::export]]
List
  V_list(double sig2_b,
          double sig2_eps,
          arma::mat ZZt_b,
          arma::mat ZZt_b_ii,
          arma::mat ZZt_eps,
          arma::mat ZZt_eps_ii,
          int m_i){
  // ZeZet = In
  arma::mat V_i = sig2_eps*ZZt_eps_ii + sig2_b*ZZt_b_ii;
  arma::mat V_inv =  inv_sympd(sig2_eps*ZZt_eps + sig2_b*ZZt_b);
  // arma::vec eval;
  // arma::mat evec;
  // eig_sym(eval, evec, V_inv(arma::span(0,m_i-1), arma::span(0,m_i-1)));

  arma::mat V_inv_1half = arma::sqrtmat_sympd(V_inv(arma::span(0,m_i-1),
                                                    arma::span(0,m_i-1)));

  // arma::mat V_inv_1half = chol(V_inv(arma::span(0,m_i-1),arma::span(0,m_i-1)));
  return List::create(Named("V_i") = V_i,
                      Named("V_inv") = V_inv,
                      Named("V_inv_1half_i") = V_inv_1half);
}


List
  V_list_lik(double sig2_b,
             double sig2_eps,
             arma::mat ZZt_b_ii,
             arma::mat ZZt_eps_ii){
    // ZeZet = In
    arma::mat V_i = sig2_eps*ZZt_eps_ii + sig2_b*ZZt_b_ii;
    arma::mat V_inv = inv_sympd(V_i);
    double ldet;
    double sign;

    arma::log_det(ldet, sign, V_i);

    return List::create(Named("V_inv") = V_inv,
                        Named("ldet") = ldet);
  }

//arma::vec repVec(arma::vec xv, int times, int xvLen){
//  arma::mat ans(times*xvLen,1);
//  ans = repmat(xv, times,1);
//  return ans;
//}

// this is REML proposal II of Richardson & Welsh 95' (RW95)
// Only two random components are allowed: random intercept and
// the error term. Psi is a vector made by (1.8) and (1.10) of RW95
//' @title Robust Restricted Likelihood Estimating Function II
//'
//' @description Robust estimating functions for a linear mixed effects model
//' based on the Restricted Likelihood. This is the robust REML proposal II
//' of Richardson & Welsh (1995).
//' @usage Psi_reml2(betas, sig2_b, sig2_eps, y, Xn, Xnt, ZZt_b,
//'             ZZt_b_ii, ZZt_eps, ZZt_eps_ii, c_hub, c2_hub, K2n,
//'             n_betas, m_i, n, n_ind)
//' @param betas The vector of fixed effects.
//' @param sig2_b The variance of the random effects.
//' @param sig2_eps The variance of the error term.
//' @param y The matrix of the response variable. Must be a ...
//' @param Xn The matrix of covariates for the whole sample
//' @param Xnt The transposed matrix of covariates for the whole sample
//' @param ZZt_b The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the random effects
//' @param ZZt_b_ii Same as \code{ZZt_b} but the block for the ith individual
//' @param ZZt_eps The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the error term.
//' @param ZZt_eps_ii Same as \code{ZZt_eps} but the block for the ith individual.
//' @param c_hub The bound for Huber's psi function for the location parameters
//' @param c2_hub The bound for Huber's proposal 2 estimating function for scale parameters.
//' @param K2n The K2 matrix consitency correction on the scale parameters.
//' @param n_betas number of fixed effects
//' @param m_i number of replications
//' @param n the sample size
//' @param n_ind the number of individuals
//' @noRd
//[[Rcpp::export]]
arma::vec
  Psi_reml2(arma::vec betas,
            double sig2_b,
            double sig2_eps,
            arma::mat y,
            arma::mat Xn,
            arma::mat Xnt,
            arma::mat ZZt_b,
            arma::mat ZZt_b_ii,
            arma::mat ZZt_eps,
            arma::mat ZZt_eps_ii,
            const double c_hub,
            const double c2_hub,
            arma::mat K2n,
            int n_betas,
            int m_i,
            int n,
            int n_ind)
  {
    // n_ind: # unique individuals
    // m_i:   # repeatitions of individual i,
    // n:     sample size
    // for balanced data only, n_ind unique individuals repeated m_i times
    // Xn is n*p x m_i design matrix equal for each individual, Xnt = t(Xn)
    // ZZt_eps & ZZt_b are n x n matrices for the random effects & errors,
    // The covariance matrix V is Sigma repeated block-wise n_ind times.

    int n_varcom = 2;
    arma::vec ans(n_betas + n_varcom, arma::fill::zeros);

    List Vstuff = V_list(sig2_b,
                         sig2_eps,
                         ZZt_b,
                         ZZt_b_ii,
                         ZZt_eps,
                         ZZt_eps_ii,
                         m_i);

    arma::mat V_inv = Vstuff["V_inv"];
    arma::mat V_inv_1half = Vstuff["V_inv_1half_i"];

    //  arma::mat Vm_in1 = inv(sig2_b*ZZt_b + sig2_eps*ZZt_eps);
    //
    // arma::vec eval;
    // arma::mat evec;
    //   eig_sym(eval, evec, Sm_in1);
    //
    //   arma::mat sqrtSm_in = evec*diagmat(sqrt(eval))*trans(evec);

    arma::vec psiHub(m_i), psiHub2(m_i), rrv(m_i);

    // matrix P for the REML
    arma::mat P = V_inv - V_inv*Xn*inv(Xnt*V_inv*Xn)*Xnt*V_inv;

    // consistency correction for the variance components
    double consistency_eps = trace(K2n*P*ZZt_eps);
    double consistency_b = trace(K2n*P*ZZt_b);

    // individual submatrix X
    arma::mat Xi(m_i, n_betas);

    for(int i=0; i<n_ind; i++){

      Xi = Xn.rows(arma::span(i*m_i, (i+1)*m_i -1));

      //vector of residuals
      rrv = V_inv_1half*(y.col(i) - Xi*betas);
      psiHub = vpsi_huber(rrv, c_hub, m_i);
      psiHub2 = vpsi_huber(rrv, c2_hub, m_i);

      //update betas
      ans(arma::span(0,n_betas-1)) += Xi.t()*V_inv_1half*psiHub;

      // update sigma_b
      ans(n_betas) += 0.5*arma::as_scalar(trans(psiHub2)*V_inv_1half*ZZt_b_ii*
        V_inv_1half*psiHub2);

      // est. eq for sigma_eps
      ans(n_betas+1) += 0.5*arma::as_scalar(trans(psiHub2)*
        V_inv(arma::span(0,m_i-1),arma::span(0,m_i-1))*psiHub2);
    }
    ans(n_betas) -= 0.5*consistency_b;
    ans(n_betas+1) -= 0.5*consistency_eps;

    //   rrv = Vn_inv_1half*(yv - Xn*beta);
    //   psiHub = vpsi_huber(rrv, c_hub, m_i*n_ind);
    //   psiHub2 = vpsi_huber(rrv, c_hub2, m_i*n_ind);
    //
    //   ans(nbeta) = 0.5*(arma::as_scalar(trans(psiHub2)*Vn_inv_1half*ZZt_b*
    //                Vn_inv_1half*psiHub2) - consistency_b);
    //   ans(nbeta+1) = 0.5*(arma::as_scalar(trans(psiHub2)*V_inv*psiHub2) - consistency_eps);
    return ans;
  }

//' @title Robust REML II Estimating Function for fixed effects
//'
//' @description Robust estimating function for the fixed effects of the linear mixed effects model
//' based on REML. This is the robust REML proposal II of Richardson & Welsh (1995).
//' @usage Psi_reml2_betas(betas, sig2_b, sig2_eps, y, Xn, ZZt_b,
//'             ZZt_b_ii, ZZt_eps, ZZt_eps_ii, c_hub,
//'             n_betas, m_i, n, n_ind)
//' @param betas The vector of fixed effects.
//' @param sig2_b The variance of the random effects.
//' @param sig2_eps The variance of the error term.
//' @param y The matrix of the response variable. Must be a ...
//' @param Xn The matrix of covariates for the whole sample
//' @param ZZt_b The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the random effects
//' @param ZZt_b_ii Same as \code{ZZt_b} but the block for the ith individual
//' @param ZZt_eps The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the error term.
//' @param ZZt_eps_ii Same as \code{ZZt_eps} but the block for the ith individual.
//' @param c_hub The bound for Huber's psi function for the location parameters
//' @param n_betas number of fixed effects
//' @param m_i number of replications
//' @param n the sample size
//' @param n_ind the number of individuals
//' @noRd
//[[Rcpp::export]]
arma::vec
  Psi_reml2_betas(arma::vec betas,
                  double sig2_b,
                  double sig2_eps,
                  arma::mat y,
                  arma::mat Xn,
                  arma::mat ZZt_b,
                  arma::mat ZZt_b_ii,
                  arma::mat ZZt_eps,
                  arma::mat ZZt_eps_ii,
                  const double c_hub,
                  int n_betas,
                  int m_i,
                  int n,
                  int n_ind)
  {
    // n_ind # unique individuals, m_i # repeatitions of individual i, n sample size
    // for balanced data only, n_ind unique individuals repeated m_i times
    // Xt is p x m_i design matrix equal for each individual, X = t(X)
    // ZZt_eps & ZZt_b are n times n matrices for the random effects & errors,
    // The covariance matrix V is Sigma repeated block-wise n_ind times.

    arma::vec ans(n_betas, arma::fill::zeros);

    List Vstuff = V_list(sig2_b,
                         sig2_eps,
                         ZZt_b,
                         ZZt_b_ii,
                         ZZt_eps,
                         ZZt_eps_ii,
                         m_i);

    arma::mat V_inv = Vstuff["V_inv"];
    arma::mat V_inv_1half = Vstuff["V_inv_1half_i"];

    arma::vec psiHub(m_i), rrv(m_i);

    arma::mat Xi(m_i, n_betas);

    for(int i=0; i<n_ind; i++){

      Xi = Xn.rows(arma::span(i*m_i, (i+1)*m_i -1));

      //vector of residuals
      rrv = V_inv_1half*(y.col(i) - Xi*betas);
      psiHub = vpsi_huber(rrv, c_hub, m_i);

      //update betas
      ans += Xi.t()*V_inv_1half*psiHub;
    }
    return ans;
  }

//' @title Robust REML II Estimating Function for the random effects
//'
//' @description Robust estimating function for the variance component of the random
//'  effects of the linear mixed effects model based on REML. This is the robust REML
//'  proposal II of Richardson & Welsh (1995).
//' @usage Psi_reml2_sig2_b(betas, sig2_b, sig2_eps, y, Xn, Xnt, ZZt_b,
//'             ZZt_b_ii, ZZt_eps, ZZt_eps_ii, c2_hub, K2n,
//'             n_betas, m_i, n, n_ind)
//' @param betas The vector of fixed effects.
//' @param sig2_b The variance of the random effects.
//' @param sig2_eps The variance of the error term.
//' @param y The matrix of the response variable. Must be a ...
//' @param Xn The matrix of covariates for the whole sample
//' @param Xnt The transposed matrix of covariates for the whole sample
//' @param ZZt_b The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the random effects
//' @param ZZt_b_ii Same as \code{ZZt_b} but the block for the ith individual
//' @param ZZt_eps The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the error term.
//' @param ZZt_eps_ii Same as \code{ZZt_eps} but the block for the ith individual.
//' @param c2_hub The bound for Huber's psi2 function for the variance paramerr of the random effects
//' @param K2n The K2 matrix consitency correction on the scale parameters.
//' @param n_betas number of fixed effects
//' @param m_i number of replications
//' @param n the sample size
//' @param n_ind the number of individuals
//' @noRd
//[[Rcpp::export]]
double
  Psi_reml2_sig2_b(arma::vec betas,
                   double sig2_b,
                   double sig2_eps,
                   arma::mat y,
                   arma::mat Xn,
                   arma::mat Xnt,
                   arma::mat ZZt_b,
                   arma::mat ZZt_b_ii,
                   arma::mat ZZt_eps,
                   arma::mat ZZt_eps_ii,
                   const double c2_hub,
                   arma::mat K2n,
                   int n_betas,
                   int m_i,
                   int n,
                   int n_ind)
  {
    // n_ind # unique individuals, m_i # repeatitions of individual i, n sample size
    // for balanced data only, n_ind unique individuals repeated m_i times
    // Xt is p x m_i design matrix equal for each individual, X = t(X)
    // ZZt_eps & ZZt_b are n times n matrices for the random effects & errors,
    // The covariance matrix V is Sigma repeated block-wise n_ind times.

    double ans = 0.0;

    List Vstuff = V_list(sig2_b,
                         sig2_eps,
                         ZZt_b,
                         ZZt_b_ii,
                         ZZt_eps,
                         ZZt_eps_ii,
                         m_i);

    arma::mat V_inv = Vstuff["V_inv"];
    arma::mat V_inv_1half = Vstuff["V_inv_1half_i"];

    //  arma::mat Vm_in1 = inv(sig2_b*ZZt_b + sig2_eps*ZZt_eps);
    //
    // arma::vec eval;
    // arma::mat evec;
    //   eig_sym(eval, evec, Sm_in1);
    //
    //   arma::mat sqrtSm_in = evec*diagmat(sqrt(eval))*trans(evec);

    arma::vec psiHub2(m_i), rrv(m_i);
    arma::mat P = V_inv - V_inv*Xn*inv(Xnt*V_inv*Xn)*Xnt*V_inv;

    double consistency_b = trace(K2n*P*ZZt_b);

    arma::mat Xi(m_i, n_betas);

    for(int i=0; i<n_ind; i++){
      Xi = Xn.rows(arma::span(i*m_i, (i+1)*m_i -1));

      //vector of residuals
      rrv = V_inv_1half*(y.col(i) - Xi*betas);
      psiHub2 = vpsi_huber(rrv, c2_hub, m_i);

      // update sigma_b
      ans += 0.5*arma::as_scalar(trans(psiHub2)*V_inv_1half*ZZt_b_ii*
        V_inv_1half*psiHub2);

    }

    return ans - 0.5*consistency_b;
  }

//' @title Robust REML II Estimating Function for the random effects
//'
//' @description Robust estimating function for the variance of the residuals
//'  of the linear mixed effects model based on REML. This is the robust REML
//'  proposal II of Richardson & Welsh (1995).
//' @usage Psi_reml2_sig2_eps(betas, sig2_b, sig2_eps, y, Xn, Xnt, ZZt_b,
//'             ZZt_b_ii, ZZt_eps, ZZt_eps_ii, c2_hub, K2n,
//'             n_betas, m_i, n, n_ind)
//' @param betas The vector of fixed effects.
//' @param sig2_b The variance of the random effects.
//' @param sig2_eps The variance of the error term.
//' @param y The matrix of the response variable. Must be a ...
//' @param Xn The matrix of covariates for the whole sample
//' @param Xnt The transposed matrix of covariates for the whole sample
//' @param ZZt_b The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the random effects
//' @param ZZt_b_ii Same as \code{ZZt_b} but the block for the ith individual
//' @param ZZt_eps The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the error term.
//' @param ZZt_eps_ii Same as \code{ZZt_eps} but the block for the ith individual.
//' @param c2_hub The bound for Huber's psi2 function for the variance paramerr of the error term
//' @param K2n The K2 matrix consitency correction on the scale parameters.
//' @param n_betas number of fixed effects
//' @param m_i number of replications
//' @param n the sample size
//' @param n_ind the number of individuals
//' @noRd
//[[Rcpp::export]]
double
  Psi_reml2_sig2_eps(arma::vec betas,
                     double sig2_b,
                     double sig2_eps,
                     arma::mat y,
                     arma::mat Xn,
                     arma::mat Xnt,
                     arma::mat ZZt_b,
                     arma::mat ZZt_b_ii,
                     arma::mat ZZt_eps,
                     arma::mat ZZt_eps_ii,
                     const double c2_hub,
                     arma::mat K2n,
                     int n_betas,
                     int m_i,
                     int n,
                     int n_ind)
  {
    // n_ind # unique individuals, m_i # repeatitions of individual i, n sample size
    // for balanced data only, n_ind unique individuals repeated m_i times
    // Xt is p x m_i design matrix equal for each individual, X = t(X)
    // ZZt_eps & ZZt_b are n times n matrices for the random effects & errors,
    // The covariance matrix V is Sigma repeated block-wise n_ind times.

    double ans = 0.0;

    List Vstuff = V_list(sig2_b,
                         sig2_eps,
                         ZZt_b,
                         ZZt_b_ii,
                         ZZt_eps,
                         ZZt_eps_ii,
                         m_i);

    arma::mat V_inv = Vstuff["V_inv"];
    arma::mat V_inv_1half = Vstuff["V_inv_1half_i"];

    //  arma::mat Vm_in1 = inv(sig2_b*ZZt_b + sig2_eps*ZZt_eps);
    //
    // arma::vec eval;
    // arma::mat evec;
    //   eig_sym(eval, evec, Sm_in1);
    //
    //   arma::mat sqrtSm_in = evec*diagmat(sqrt(eval))*trans(evec);

    arma::vec psiHub2(m_i), rrv(m_i);
    arma::mat P = V_inv - V_inv*Xn*inv(Xnt*V_inv*Xn)*Xnt*V_inv;

    double consistency_eps = trace(K2n*P*ZZt_eps);

    arma::mat Xi(m_i, n_betas);

    for(int i=0; i<n_ind; i++){

      Xi = Xn.rows(arma::span(i*m_i, (i+1)*m_i -1));
      //vector of residuals
      rrv = V_inv_1half*(y.col(i) - Xi*betas);
      psiHub2 = vpsi_huber(rrv, c2_hub, m_i);

      // update sigma_eps
      ans += 0.5*arma::as_scalar(trans(psiHub2)*
            V_inv(arma::span(0,m_i-1),arma::span(0,m_i-1))*psiHub2);

    }

    return ans - 0.5*consistency_eps;
  }

//same as as _reml2 but without computing the (expensive)
//consistency corrections

//[[Rcpp::export]]
arma::vec
  Psi_rlmm_reml2_abc(arma::vec betas,
                     arma::mat V_inv,
                     arma::mat V_inv_1half,
                     arma::mat y,
                     arma::mat Xn,
                     arma::mat Xnt,
                     arma::mat ZZt_b_ii,
                     const double c_hub,
                     const double c2_hub,
                     int n_betas,
                     int m_i,
                     int n,
                     int n_ind,
                     double consistency_b,
                     double consistency_eps)
  {
    // n_ind # unique individuals, m_i # repeatitions of individual i, n sample size
    // for balanced data only, n_ind unique individuals repeated m_i times
    // Xt is p x m_i design matrix equal for each individual, X = t(X)
    // ZZt_eps & ZZt_b are n times n matrices for the random effects & errors,
    // The covariance matrix V is Sigma repeated block-wise n_ind times.

    arma::vec ans(n_betas+2, arma::fill::zeros);

//     List Vstuff = V_list_abc(sig2_b, sig2_eps, ZZt_b_ii, ZZt_eps_ii);
//
//     arma::mat V_inv = Vstuff["V_inv"];
//     arma::mat V_inv_1half = Vstuff["V_inv_1half_i"];

    arma::vec psiHub(m_i), psiHub2(m_i), rrv(m_i);
    arma::mat Xi(m_i, n_betas);

    for(int i=0; i<n_ind; i++){
      Xi = Xn.rows(arma::span(i*m_i, (i+1)*m_i -1));
      //vector of residuals
      rrv = V_inv_1half*(y.col(i) - Xi*betas);
      psiHub = vpsi_huber(rrv, c_hub, m_i);
      psiHub2 = vpsi_huber(rrv, c2_hub, m_i);

      //update betas
      ans(arma::span(0,n_betas-1)) += Xi.t()*V_inv_1half*psiHub;

      // update sigma_b
      ans(n_betas) += 0.5*arma::as_scalar(trans(psiHub2)*V_inv_1half*ZZt_b_ii*
        V_inv_1half*psiHub2);

      // est. eq for sigma_eps
      ans(n_betas+1) += 0.5*arma::as_scalar(trans(psiHub2)*V_inv*psiHub2);
    }

    ans(n_betas) -= 0.5*consistency_b;
    ans(n_betas+1) -= 0.5*consistency_eps;

    return ans;
  }

//[[Rcpp::export]]
arma::vec
  Psi_rlmm_reml2_ith(arma::vec betas,
                  double sig2_b,
                  double sig2_eps,
                  arma::mat y,
                  arma::mat Xn,
                  arma::mat Xnt,
                  arma::mat ZZt_b,
                  arma::mat ZZt_b_ii,
                  arma::mat ZZt_eps,
                  arma::mat ZZt_eps_ii,
                  const double c_hub,
                  const double c2_hub,
                  arma::mat K2n,
                  int n_betas,
                  int m_i,
                  int n,
                  int n_ind,
                  int ith)
  {
    // n_ind # unique individuals, m_i # repeatitions of individual i, n sample size
    // for balanced data only, n_ind unique individuals repeated m_i times
    // Xt is p x m_i design matrix equal for each individual, X = t(X)
    // ZZt_eps & ZZt_b are n times n matrices for the random effects & errors,
    // The covariance matrix V is Sigma repeated block-wise n_ind times.

    arma::vec ans(n_betas+2);

    List Vstuff = V_list(sig2_b,
                         sig2_eps,
                         ZZt_b,
                         ZZt_b_ii,
                         ZZt_eps,
                         ZZt_eps_ii,
                         m_i);

    arma::mat V_inv = Vstuff["V_inv"];
    arma::mat V_inv_1half = Vstuff["V_inv_1half_i"];

    arma::vec psiHub(m_i), psiHub2(m_i), rrv(m_i);
    arma::mat P = V_inv - V_inv*Xn*inv(Xnt*V_inv*Xn)*Xnt*V_inv;

    double consistency_eps = accu(diagmat(K2n*P*ZZt_eps));
    double consistency_b = accu(diagmat(K2n*P*ZZt_b));

    // for(int i=0; i<n_ind; i++){
      //vector of residuals
    arma::mat Xith = Xn.rows(arma::span(ith*m_i, (ith+1)*m_i -1));
    rrv = V_inv_1half*(y.col(ith-1) - Xith*betas);
    psiHub = vpsi_huber(rrv, c_hub, m_i);
    psiHub2 = vpsi_huber(rrv, c2_hub, m_i);

      //update betas
    ans(arma::span(0,n_betas-1)) = Xith.t()*V_inv_1half*psiHub;

      // update sigma_b
    ans(n_betas) = 0.5*(arma::as_scalar(arma::trans(psiHub2)*V_inv_1half*ZZt_b_ii*
        V_inv_1half*psiHub2) - consistency_b/n_ind);

      // est. eq for sigma_eps
    ans(n_betas+1) = 0.5*(arma::as_scalar(trans(psiHub2)*
                        V_inv(arma::span(0,m_i-1),arma::span(0,m_i-1))*psiHub2) -
                        consistency_eps/n_ind);
  return ans;
}

//[[Rcpp::export]]
arma::mat
  Psi_rlmm_reml2_alli(arma::vec betas,
                     double sig2_b,
                     double sig2_eps,
                     arma::mat y,
                     arma::mat Xn,
                     arma::mat Xnt,
                     arma::mat ZZt_b,
                     arma::mat ZZt_b_ii,
                     arma::mat ZZt_eps,
                     arma::mat ZZt_eps_ii,
                     const double c_hub,
                     const double c2_hub,
                     arma::mat K2n,
                     int n_betas,
                     int m_i,
                     int n,
                     int n_ind)
  {
    // n_ind # unique individuals, m_i # repeatitions of individual i, n sample size
    // for balanced data only, n_ind unique individuals repeated m_i times
    // Xt is p x m_i design matrix equal for each individual, X = t(X)
    // ZZt_eps & ZZt_b are n times n matrices for the random effects & errors,
    // The covariance matrix V is Sigma repeated block-wise n_ind times.

    arma::mat ans(n_betas+2, n_ind);

    List Vstuff = V_list(sig2_b,
                         sig2_eps,
                         ZZt_b,
                         ZZt_b_ii,
                         ZZt_eps,
                         ZZt_eps_ii,
                         m_i);

    arma::mat V_inv = Vstuff["V_inv"];
    arma::mat V_inv_1half = Vstuff["V_inv_1half_i"];

    arma::vec psiHub(m_i), psiHub2(m_i), rrv(m_i);
    arma::mat P = V_inv - V_inv*Xn*inv(Xnt*V_inv*Xn)*Xnt*V_inv;

    double consistency_eps = trace(K2n*P*ZZt_eps);
    double consistency_b = trace(K2n*P*ZZt_b);

    arma::mat Xi(m_i, n_betas);

    for(int i=0; i<n_ind; i++){

      Xi = Xn.rows(arma::span(i*m_i, (i+1)*m_i -1));
    //vector of residuals
      rrv = V_inv_1half*(y.col(i) - Xi*betas);
      psiHub = vpsi_huber(rrv, c_hub, m_i);
      psiHub2 = vpsi_huber(rrv, c2_hub, m_i);

    //update betas
      ans(arma::span(0,n_betas-1), i) += Xi.t()*V_inv_1half*psiHub;

    // update sigma_b
      ans(n_betas, i) = 0.5*(arma::as_scalar(arma::trans(psiHub2)*V_inv_1half*ZZt_b_ii*
      V_inv_1half*psiHub2) - consistency_b/n_ind);

    // est. eq for sigma_eps
      ans(n_betas+1, i) = 0.5*(arma::as_scalar(trans(psiHub2)*
                        V_inv(arma::span(0,m_i-1),arma::span(0,m_i-1))*psiHub2) -
                        consistency_eps/n_ind);
    }
  return ans;
}


// this is (1.7) of Richardson & Welsh 95'
// arma::vec funPsi_rlmm(arma::vec beta,
//                 double sig2_b,
//                 double sig2_eps,
//                 arma::mat ymat,
//                 arma::mat X,
//                 arma::mat Xt,
//                 arma::mat ZZt_b,
//                 arma::mat ZZtj_b,
//                 arma::mat ZZt_eps,
//                 arma::mat Xn,
//                 arma::mat Xnt,
//                 const double c_hub,
//                 const double c_hub2,
//                 arma::mat K2,
//                 int nbeta,
//                 int m_i,
//                 int n,
//                 int n_ind)
// {
//   // n_ind # unique individuals, m_i # repeatitions of individual i, n sample size
//   // data are balanced, n_ind unique individuals repeated m_i times
//   // Xt is p x m_i design matrix equal for each individual, X = t(X)
//   // ZZt_eps & ZZt_b are n times n matrices for the random effects & errors,
//   // The covariance matrix V is Sigma repeated n_ind times.
//
//   arma::vec ans(nbeta+2, arma::fill::zeros);
//
//   arma::mat Sm_in1(m_i, m_i);
//   Sm_in1 = Smat_inv(sig2_b, sig2_eps, Sm_in1);
//
// //  arma::mat Vm_in1 = inv(sig2_b*ZZt_b + sig2_eps*ZZt_eps);
// //
//   arma::vec eval;
//   arma::mat evec;
//   eig_sym(eval, evec, Sm_in1);
//
//   arma::mat sqrtSm_in = evec*diagmat(sqrt(eval))*trans(evec);
//   arma::vec psiHub(m_i);
//   arma::vec psiHub2(m_i);
//   arma::vec rrv(m_i);
// //  arma::mat P = Vm_in1 - Vm_in1*Xn*inv(Xnt*Vm_in1*Xn)*Xnt*Vm_in1;
//
//   double ccor_eps = accu(diagmat(K2*Sm_in1));
//   double ccor_b = accu(diagmat(K2*Sm_in1*ZZtj_b));
//
//   for(int i=0; i<n_ind; i++){
//     //vector of residuals
//     rrv = sqrtSm_in*(trans(ymat.row(i)) - X*beta);
//     psiHub = vpsi_hub(rrv, c_hub, m_i);
//     psiHub2 = vpsi_hub(rrv, c_hub2, m_i);
//     //update betas
//     ans(arma::span(0,nbeta-1)) += Xt*sqrtSm_in*psiHub;
//
//     // update sigma_b
//     ans(nbeta) += 0.5*(as_scalar(trans(psiHub2)*sqrtSm_in*ZZtj_b*
//                     sqrtSm_in*psiHub2) - ccor_b);
//
//   // est. eq for sigma_eps
//     ans(nbeta+1) += 0.5*(as_scalar(trans(psiHub2)*Sm_in1*psiHub2) - ccor_eps);
//   }
//   return ans;
// }

// //[[Rcpp::export]]
// arma::mat funPsi_rlmm_aux(arma::vec beta,
//                 double sig2_b,
//                 double sig2_eps,
//                 arma::mat ymat,
//                 arma::mat X,
//                 arma::mat Xt,
//                 arma::mat ZZt_b,
//                 arma::mat ZZtj_b,
//                 arma::mat ZZt_eps,
//                 arma::mat Xn,
//                 arma::mat Xnt,
//                 const double c_hub,
//                 const double c_hub2,
//                 arma::mat K2,
//                 int nbeta,
//                 int m_i,
//                 int n,
//                 int n_ind)
// {
//   // n_ind # unique individuals, m_i # repeatitions of individual i, n sample size
//   // data are balanced, n_ind unique individuals repeated m_i times
//   // Xt is p x m_i design matrix equal for each individual, X = t(X)
//   // ZZt_eps & ZZt_b are n times n matrices for the random effects & errors,
//   // The covariance matrix V is Sigma repeated n_ind times.
//
//   arma::mat ans(nbeta+2, n_ind, arma::fill::zeros);
//   arma::mat Sm_in1(m_i, m_i);
//   Sm_in1 = Smat_inv(sig2_b, sig2_eps, Sm_in1);
// //  arma::mat Vm_in1 = inv(sig2_b*ZZt_b + sig2_eps*ZZt_eps);
// //  arma::mat cholVm_in1 = trans(chol(Vm_in1));
//   arma::vec eval;
//   arma::mat evec;
//   eig_sym(eval, evec, Sm_in1);
//   arma::mat sqrtSm_in = evec*diagmat(sqrt(eval))*trans(evec);
//   arma::vec psiHub(m_i);
//   arma::vec psiHub2(m_i);
//   arma::vec rrv(m_i);
// //  arma::mat P = Vm_in1 - Vm_in1*Xn*inv(Xnt*Vm_in1*Xn)*Xnt*Vm_in1;
//
// //  double ccor_eps = accu(diagmat(K2*P*ZZt_eps));
// //  double ccor_b = accu(diagmat(K2*P*ZZt_b));
//
//   double ccor_eps = accu(diagmat(K2*Sm_in1));
//   double ccor_b = accu(diagmat(K2*Sm_in1*ZZtj_b));
// //
//   for(int i=0; i<n_ind; i++){
//     //vector of residuals
// //    psiv = vpsi_biwgt(sqrtSm_in*(trans(ymat.row(i)) - X*beta), c_tune, m_i);
//     rrv = sqrtSm_in*(trans(ymat.row(i)) - X*beta);
//     psiHub = vpsi_hub(rrv, c_hub, m_i);
//     psiHub2 = vpsi_hub(rrv, c_hub2, m_i);
//
//     //update betas
//     ans(arma::span(0,nbeta-1),i) = Xt*sqrtSm_in*psiHub;
//
//     // update sigma_b
//     ans(nbeta,i) = 0.5*(as_scalar(trans(psiHub2)*sqrtSm_in*ZZtj_b*
//                     sqrtSm_in*psiHub2) - ccor_b);
//
//   // est. eq for sigma_eps
//   ans(nbeta+1, i) = 0.5*(as_scalar(trans(psiHub2)*Sm_in1*psiHub2) - ccor_eps);
//   }
//   return ans;
// }

//[[Rcpp::export]]
arma::mat
simData_rlmm(arma::vec betas,
             double sig2_b,
             double sig2_eps,
             arma::mat Xn,
             arma::mat ZZt_b_ii,
             arma::mat ZZt_eps_ii,
             int m_i,
             int n_ind)
{
 arma::mat ans(m_i, n_ind),
 V_i = sig2_eps*ZZt_eps_ii + sig2_b*ZZt_b_ii;
 arma::mat A;
 chol(A, V_i); // note: this gives A^TA = S, with A upper-triangular
 arma::mat At = A.t();

 bool status = A.is_empty();

 if(status == true){

   Rprintf("\nCholesky decomposition in simData_rlmm failed to converge!");
   return arma::mat(arma::datum::nan, arma::datum::nan);

 } else {

   RNGScope scope;
   arma::vec mui(m_i);

   for(int i=0; i<n_ind; i++){

     mui = Xn.rows(arma::span(i*m_i, (i+1)*m_i -1))*betas;
     ans.col(i) = rmvnorm2(mui, At, m_i);
     }
   return ans;
   }
}

// //[[Rcpp::export]]
// arma::mat
// ABC_rlmm_reml2(int nabc,
//                double eps,
//                arma::mat y,
//                arma::mat X,
//                arma::mat ZZt_b,
//                arma::mat ZZt_b_ii,
//                arma::mat ZZt_eps,
//                arma::mat ZZt_eps_ii,
//                arma::mat Xn,
//                const double c_hub,
//                const double c2_hub,
//                arma::mat K2n,
//                arma::vec param_hat,
//                arma::mat prop_scale_mat,
//                arma::vec Psi_hat,
//                arma::mat J_Psi_hat,
//                double prior_beta_stdev,
//                double prior_sig2_scale,
//                arma::vec thin,
//                int trace_int)
// {
//   int nacc = 1,
//   n_param = param_hat.n_elem,
//   n_betas = n_param - 2,
//   m_i = X.n_rows, n = Xn.n_rows, n_ind = y.n_cols;
//   double sig2_b_hat = exp(param_hat(n_betas)),
//   sig2_eps_hat = exp(param_hat(n_betas+1));
//
//   arma::mat Xt = trans(X),
//   Xnt = trans(Xn),
//   ans(n_param, nabc),
//   y_sim(m_i, n_ind);
//
//   arma::vec param_prop(n_param),
//   Psi_sim(n_param),
//   betas_hat = param_hat(arma::span(0, n_betas-1));
//
//
//   List Vstuff_hat = V_list(sig2_b_hat, sig2_eps_hat, ZZt_b, ZZt_b_ii,
//                            ZZt_eps, ZZt_eps_ii, m_i);
//
//   arma::mat V_inv_hat = Vstuff_hat["V_inv"],
//   V_inv_1half_hat = Vstuff_hat["V_inv_1half_i"],
//   P_hat = V_inv_hat - V_inv_hat*Xn*inv(Xnt*V_inv_hat*Xn)*Xnt*V_inv_hat,
//   V_inv_hat_ii = V_inv_hat(arma::span(0,m_i-1), arma::span(0,m_i-1));
//
//   double consistency_eps_hat = trace(K2n*P_hat*ZZt_eps),
//   consistency_b_hat = trace(K2n*P_hat*ZZt_b);
//
//   double prop_df = 5.0, rho_dist = 0.0, log_accept_ratio = 0.0;
//
// /*start algorithm from hat.theta*/
//   ans.col(0) = param_hat;
//
// /* for i to nabc do*/
//   RNGScope scope;
//   for(int i = 1; i < nabc; i++){
//     // draw from the jumping density
//     param_prop = rmvt(ans.col(i-1), prop_scale_mat, n_param, prop_df);
//
//     // acceptance ratio of prior and jumping density
//     log_accept_ratio = dPrior_lmm(param_prop(arma::span(0,n_betas-1)),
//                                   param_prop(n_betas),
//                                   param_prop(n_betas+1),
//                                   prior_beta_stdev, prior_sig2_scale, n_betas, true) -
//                        dPrior_lmm(ans(arma::span(0,n_betas-1), i-1),
//                                   ans(n_betas, i-1),
//                                   ans( n_betas+1, i-1),
//                                   prior_beta_stdev, prior_sig2_scale, n_betas, true);
//
//
//     if (log(as_scalar(arma::randu(1))) > log_accept_ratio) {
//       ans.col(i) = ans.col(i-1);
//       } else {
//         // compute the ABC disntace
//     y_sim = simData_rlmm(param_prop(arma::span(0,n_betas-1)),
//                             exp(param_prop(n_betas)),
//                             exp(param_prop(n_betas+1)),
//                             X, ZZt_b_ii, ZZt_eps_ii,
//                             m_i, n_ind);
// //     Psi_sim = Psi_rlmm_reml2(betas_hat, sig2_b_hat, sig2_eps_hat,
// //                                  y_sim, X, Xt, ZZt_b,
// //                                  ZZt_b_ii, ZZt_eps,
// //                                  ZZt_eps_ii, Xn, Xnt,
// //                                  c_hub, c2_hub, K2n = K2n,
// //                                  n_betas, m_i, n, n_ind);
//
//     Psi_sim = Psi_rlmm_reml2_abc(betas_hat, V_inv_hat_ii, V_inv_1half_hat,
//                                  y_sim, X, Xt, ZZt_b_ii, c_hub, c2_hub,
//                                  n_betas, m_i, n, n_ind,
//                                  consistency_b_hat, consistency_eps_hat);
//
//     rho_dist = arma::as_scalar(trans(Psi_sim-Psi_hat)*J_Psi_hat*(Psi_sim-Psi_hat));
//     if (rho_dist > eps){
//        // take the previous value
//       ans.col(i) = ans.col(i-1);
//     } else {
//       // accept
//       nacc += 1;
//       ans.col(i) = param_prop;
//     }
//   }
//   // print some information on the simulation process
//   if(i%trace_int == 0 ) Rprintf("\rABC %d of %d, accpted %d, acc.ratio %f\r", i,
//      nabc, nacc, static_cast<double>(nacc)/static_cast<double>(i));
//
//   // for stopping R computations
//   R_CheckUserInterrupt();
//   }
//
//   Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
//   Rprintf("ABC-MCMC accepteance ratio was: %3.5f",
//         static_cast<double>(nacc) / static_cast<double>(nabc));
//   Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
//
//   return thinMat(ans, thin);
//   // return rho_distv;
// }

//[[Rcpp::export]]
List
  ABCkern_reml2(int nabc,
                     double h,
                     arma::mat y,
                     arma::mat Xn,
                     arma::mat ZZt_b,
                     arma::mat ZZt_b_ii,
                     arma::mat ZZt_eps,
                     arma::mat ZZt_eps_ii,
                     const double c_hub,
                     const double c2_hub,
                     arma::mat K2n,
                     arma::vec param_hat,
                     arma::mat lch_prop_scale_mat,
                     arma::vec Psi_hat,
                     arma::mat chol_J_Psi_hat,
                     double prior_beta_stdev,
                     double prior_sig2_scale,
                     arma::vec thin,
                     int trace_int,
                     int m_i,
                     int n_ind)
  {
    int nacc = 1,
      n_param = param_hat.n_elem,
      n_betas = n_param - 2,
      n = Xn.n_rows;

    double sig2_b_hat = exp(param_hat(n_betas)),
      sig2_eps_hat = exp(param_hat(n_betas+1));

    arma::mat Xnt = trans(Xn),
      ans(n_param, nabc),
      y_prop(m_i, n_ind);

    arma::vec param_prop(n_param),
      Psi_prop(n_param),
      betas_hat = param_hat(arma::span(0, n_betas-1)),
      log_kern_rhov(nabc), scaled_Psi(n_param);


    List Vstuff_hat = V_list(sig2_b_hat,
                             sig2_eps_hat,
                             ZZt_b,
                             ZZt_b_ii,
                             ZZt_eps,
                             ZZt_eps_ii,
                             m_i);

    arma::mat V_inv_hat = Vstuff_hat["V_inv"],
      V_inv_1half_hat = Vstuff_hat["V_inv_1half_i"],
      P_hat = V_inv_hat - V_inv_hat*Xn*inv(Xnt*V_inv_hat*Xn)*Xnt*V_inv_hat,
      V_inv_hat_ii = V_inv_hat(arma::span(0,m_i-1), arma::span(0,m_i-1));

    double consistency_eps_hat = trace(K2n*P_hat*ZZt_eps),
      consistency_b_hat = trace(K2n*P_hat*ZZt_b),
      prop_df = 5.0, log_accept_ratio = 0.0, log_kern_rho_ratio = 0.0,
      log_prior_ratio = 0.0, log_kern_rho_prop = 0.0;

    /*start algorithm from hat.theta*/
    ans.col(0) = param_hat;
    log_kern_rhov(0) = n_param*Rf_dnorm4(0.0, 0.0, 1.0, true);
    // log_kern_rhov(0) = dmvnorm_I(0.0, h, n_param, true)

    /* for i to nabc do*/

    RNGScope scope;
    for(int i = 1; i < nabc; i++){
      // draw from the jumping density
      param_prop = rmvt(ans.col(i-1), lch_prop_scale_mat, n_param, prop_df);

      // acceptance ratio of prior and jumping density
      log_prior_ratio = dPrior_lmm(param_prop(arma::span(0,n_betas-1)),
                                    param_prop(n_betas),
                                    param_prop(n_betas+1),
                                    prior_beta_stdev, prior_sig2_scale, n_betas, true) -
                         dPrior_lmm(ans(arma::span(0,n_betas-1), i-1),
                                    ans(n_betas, i-1),
                                    ans(n_betas+1, i-1),
                                    prior_beta_stdev, prior_sig2_scale, n_betas, true);

        // draw data from the
      y_prop = simData_rlmm(param_prop(arma::span(0,n_betas-1)),
                           exp(param_prop(n_betas)),
                           exp(param_prop(n_betas+1)),
                           Xn,
                           ZZt_b_ii,
                           ZZt_eps_ii,
                           m_i,
                           n_ind);
      if(y_prop.is_empty() == true){
        i = i - 1;
        } else {
        //     Psi_sim = Psi_rlmm_reml2(betas_hat, sig2_b_hat, sig2_eps_hat,
        //                                  y_sim, X, Xt, ZZt_b,
        //                                  ZZt_b_ii, ZZt_eps,
        //                                  ZZt_eps_ii, Xn, Xnt,
        //                                  c_hub, c2_hub, K2n = K2n,
        //                                  n_betas, m_i, n, n_ind);
        // compute the summary statistic
      Psi_prop = Psi_rlmm_reml2_abc(betas_hat,
                                    V_inv_hat_ii,
                                    V_inv_1half_hat,
                                    y_prop,
                                    Xn,
                                    Xnt,
                                    ZZt_b_ii,
                                    c_hub,
                                    c2_hub,
                                    n_betas,
                                    m_i,
                                    n,
                                    n_ind,
                                    consistency_b_hat,
                                    consistency_eps_hat);
        // compute the distance
      scaled_Psi = chol_J_Psi_hat*(Psi_prop-Psi_hat);
        // compute the ratio of kernels
      log_kern_rho_prop = dmvnorm_I(scaled_Psi, h, n_param, true);
      log_kern_rho_ratio = log_kern_rho_prop - log_kern_rhov(i-1);
      // compute the log accetance ratio

      log_accept_ratio = log_prior_ratio + log_kern_rho_ratio;
        // pass_to_min(1) = exp(log_accept_ratio);

        // accept_ratio = min(pass_to_min);

      ans.col(i) = ans.col(i-1);
      log_kern_rhov(i) = log_kern_rhov(i-1);

      if (log(Rf_runif(0.0, 1.0)) <= log_accept_ratio){
          nacc += 1;
          ans.col(i) = param_prop;
          log_kern_rhov(i) = log_kern_rho_prop;
      }

      // print some information on the simulation process
      if(i%trace_int == 0 ) Rprintf("\rABC %d of %d, accpted %d, acc.ratio %f\r", i,
         nabc, nacc, static_cast<double>(nacc)/static_cast<double>(i));

      // for stopping R computations
      R_CheckUserInterrupt();
    }
    }
    // print final acceptance ratio
//     Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
//     Rprintf("ABC-MCMC accepteance ratio was: %3.5f",
//             static_cast<double>(nacc) / static_cast<double>(nabc));
//     Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
//     int thin = round(nabc/nacc);
//     arma::vec thinv =
    return List::create(Named("abc") = thinMat(ans, thin),
                        Named("effi") =  static_cast<double>(nacc)/
                          static_cast<double>(nabc));
    // return rho_distv;
  }

//[[Rcpp::export]]
double
log_lik_lmm(arma::vec betas,
          double sig2_b,
          double sig2_eps,
          arma::mat y,
          arma::mat Xn,
          arma::mat ZZt_b_ii,
          arma::mat ZZt_eps_ii,
          int n_betas,
          int m_i,
          int n,
          int n_ind){
  List Vstuff = V_list_lik(sig2_b,
                           sig2_eps,
                           ZZt_b_ii,
                           ZZt_eps_ii);
  double ldet = Vstuff["ldet"];
  arma::mat V_inv = Vstuff["V_inv"];
  double ans = 0.0;

  arma::mat Xi(m_i, n_betas);

  for(int i=0; i<n_ind; i++){
    Xi = Xn.rows(arma::span(i*m_i, (i+1)*m_i -1));
    ans += ldet + arma::as_scalar(trans(y.col(i)-
      Xi*betas)*V_inv*(y.col(i)-Xi*betas));
  }
  return -0.5*ans;
}
