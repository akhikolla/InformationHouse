/*
 * This code computes the asymptotic variances of these VUS estimates.
 * We use the RcppArmadillo package to compute the transpose and inverse of
 * matrix or vector.
 */

#define ARMA_NO_DEBUG
#define ARMA_USE_BLAS

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// indicate function
static inline double indvus(double a, double b, double c) {
  if((a < b) && (b < c)) {
    return 1.0;
  }else if((a < b) && (b == c)){
    return 0.5;
  }else if((a == b) && (b < c)){
    return 0.5;
  }else if((a == b) && (b == c)){
    return 1.0/6;
  }
  else{
    return 0.0;
  }
}

//==============================================================================//
/*
 * estimating function of VUS
 * Note:
 * dd_1i, dd_2j, dd_3k are estimated verison of D_1i, D_2j and D_3k, respectively.
 * E.g., FI: dd_1i, dd_2j and dd_3k correspond to rho_1i, rho_2j and rho_3k.
 * bias_ijk = I_ijk - mu_hat
 */

static inline double vusEstFunc(double dd_1i, double dd_2j, double dd_3k,
                                double bias_ijk){
  return dd_1i*dd_2j*dd_3k*bias_ijk;
}

//==============================================================================//
/*
 * derivative of estimating function of VUS
 * Note:
 * deri_dd_1i, deri_dd_2j and deri_dd_3k are the derivatives of estimated verision
 * of D_1i, D_2j and D_3k.
 */

static inline vec vusDerEstFunc(double dd_1i, double dd_2j, double dd_3k,
                                vec deri_dd_1i, vec deri_dd_2j, vec deri_dd_3k,
                                double bias_ijk){
  vec tem1 = deri_dd_1i*dd_2j*dd_3k;
  tem1 += dd_1i*deri_dd_2j*dd_3k;
  tem1 += dd_1i*dd_2j*deri_dd_3k;
  return tem1*bias_ijk;
}

//==============================================================================//
// Asymptotic Variance of VUS Estimators

//[[Rcpp::export]]
NumericVector asyVarVUS_C(NumericVector tt, NumericMatrix D_hat, double mu_hat,
                          NumericMatrix EstFunc, NumericMatrix Hess_inv,
                          NumericMatrix Der_D1_hat, NumericMatrix Der_D2_hat,
                          NumericMatrix Der_D3_hat){
  vec T = as<vec>(tt);
  mat Dhat = as<mat>(D_hat);
  mat SS = as<mat>(Hess_inv);
  mat U2 = as<mat>(EstFunc);
  mat der_Dhat1 = as<mat>(Der_D1_hat);
  mat der_Dhat2 = as<mat>(Der_D2_hat);
  mat der_Dhat3 = as<mat>(Der_D3_hat);
  uword nn = T.n_elem;
  vec S1(nn, fill::zeros), U1(SS.n_rows, fill::zeros);
  for(uword i = 0; i < nn; ++i){
    double tmp2 = 0.0;
    double t_i = T[i];
    rowvec Dhat_i = Dhat.row(i);
    for(uword j = 0; j < nn; ++j){
      if(j != i){
        double t_j = T[j];
        rowvec Dhat_j = Dhat.row(j);
        for(uword k = 0; k < nn; ++k){
          if((k != j) && (k != i)){
            double t_k = T[k];
            double i_ijk = indvus(t_i, t_j, t_k) - mu_hat;
            double i_jik = indvus(t_j, t_i, t_k) - mu_hat;
            double i_kji = indvus(t_k, t_j, t_i) - mu_hat;
            rowvec Dhat_k = Dhat.row(k);
            double tmp1 = vusEstFunc(Dhat_i[0], Dhat_j[1], Dhat_k[2], i_ijk);
            tmp1 += vusEstFunc(Dhat_j[0], Dhat_i[1], Dhat_k[2], i_jik);
            tmp1 += vusEstFunc(Dhat_k[0], Dhat_j[1], Dhat_i[2], i_kji);
            tmp2 += tmp1;
            U1 += vusDerEstFunc(Dhat_i[0], Dhat_j[1], Dhat_k[2], der_Dhat1.col(i),
                                der_Dhat2.col(j), der_Dhat3.col(k), i_ijk);
          }
        }
      }
    }
    S1[i] = tmp2;
  }
  vec Q_vus = (S1 + U2*SS*U1)/((nn - 1)*(nn - 2));
  return wrap(Q_vus);
}
// END
