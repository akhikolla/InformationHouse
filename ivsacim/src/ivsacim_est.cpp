# include <Rcpp.h>
using namespace Rcpp;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' This is the main function to compute the estimator in C++
//' @keywords internal
//' @export
// [[Rcpp::export]]
SEXP ivsacim_est(NumericVector time,
                 NumericVector event,
                 NumericVector stime, 
                 NumericVector Zc, 
                 NumericMatrix D_status,
                 NumericMatrix eps_2,
                 NumericMatrix Zc_dot){
  
  // variables for estimation
  int n = time.size();
  int k = stime.size();
  NumericVector B_D(k), dB_D(k), b_numer_D(k), b_dinom_D(k), risk_cumsum(k), indik_v(k);
  NumericMatrix dN(n, k),  risk_t(n, k), B_intD(n, k + 1);
  double tmp_dinom_D;
  double beta = 0, tot_risk = 0;
  double del = 0.01;
  
  // variables for variance estimate
  int pdim = Zc_dot.ncol();
  NumericVector B_D_se(k);
  //double B_D_se_tmp = 0;
  NumericMatrix b_var_est(k, k);
  NumericMatrix eps_1(n, k), eps(n, k);
  NumericMatrix tmp_eps_1(n, k), tmp_eps_2(n, k);
  NumericMatrix b_D_dot(pdim, k);
  NumericMatrix H_dot(k, k);
  NumericMatrix H_numer_dot(k, k), H_denom_dot(k, k), H_dot_Z(n, k);
  NumericVector eps_beta(n);
  double tmp, beta_se = 0;
  
  for (int j = 0; j < k; j++) {
    
    for (int i = 0; i < n; i++) {
      // estimation
      dN(i, j) = (time[i] == stime[j]) ? 1:0;
      dN(i, j) *= event[i];
      risk_t(i, j) = (time[i] >= stime[j]) ? 1:0;
      risk_cumsum[j] += risk_t(i, j);
      b_numer_D[j] += Zc[i] * exp(B_intD(i, j)) * dN(i, j);
      b_dinom_D[j] += Zc[i] * risk_t(i, j) * exp(B_intD(i, j)) * D_status(i, j);
    }
    
    
    tmp_dinom_D = b_dinom_D[j];
    tmp_dinom_D = (tmp_dinom_D < 0) ? -tmp_dinom_D:tmp_dinom_D;
    indik_v[j] = (tmp_dinom_D < del) ? 0:1;
    if (indik_v[j]) {
      dB_D[j] = b_numer_D[j] / b_dinom_D[j];
    }
    
    
    // estimation
    for (int i = 0; i < n; i++) {
      B_intD(i, j + 1) += D_status(i, j) * dB_D[j];
      B_intD(i, j + 1) += B_intD(i, j);
      eps_1(i, j) += Zc[i] * exp(B_intD(i, j)) * (dN(i, j) - risk_t(i, j) * D_status(i, j) * dB_D[j]);
    }
    
    
    // variance estimate
    if (indik_v[j]) {
      for (int i = 0; i < n; i++) {
        eps_1(i, j) /= b_dinom_D[j];
      }
    }
    else{
      for (int i = 0; i < n; i++) {
        eps_1(i, j) = 0;
      }
    }
    
    
    // estimation part
    beta += risk_cumsum[j] * dB_D[j];
    
    if (j == 0) {
      B_D[j] = 0 + dB_D[j];
      tot_risk += risk_cumsum[j] * (stime[j] - 0);
    }
    else{
      B_D[j] = B_D[j - 1] + dB_D[j];
      tot_risk += risk_cumsum[j] * (stime[j] - stime[j - 1]);
    }
    
  }
  
  // time invariant intensity
  beta /= tot_risk;
  
  // variance estimate
  for (int j = 0; j < k; j++) {
    for (int l = 0; l < j; l++) {
      for (int i = 0; i < n; i++) {
        H_numer_dot(l, j) += Zc[i] * exp(B_intD(i, j)) * dN(i, j) * D_status(i, l);
        H_denom_dot(l, j) += Zc[i] * risk_t(i, j) * exp(B_intD(i, j)) * D_status(i, j) * D_status(i, l);
      }
      if (indik_v[j]) {
        H_dot(l, j) = H_numer_dot(l, j) / b_dinom_D[j] - H_denom_dot(l, j) * b_numer_D[j] / b_dinom_D[j] / b_dinom_D[j];
      }
    }
  }
  
  
  
  for (int j = 0; j < k; j++) {
    for (int j1 = 0; j1 < pdim; j1++){
      for (int l = 0; l < j; l++) {
        b_D_dot(j1, j) += H_dot(l, j) * b_D_dot(j1, l);
      }
      for (int i = 0; i < n; i++) {
        if (indik_v[j]) {
          H_dot_Z(i, j) += exp(B_intD(i, j)) * dN(i, j) / b_dinom_D[j];
          H_dot_Z(i, j) -= risk_t(i, j) * exp(B_intD(i, j)) * D_status(i, j) * b_numer_D[j] / b_dinom_D[j] / b_dinom_D[j];
        }
        b_D_dot(j1, j) -= H_dot_Z(i, j) * Zc_dot(i, j1);
      }
    }
  }
  
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      for (int l = 0; l < j; l++) {
        tmp -= H_dot(l, j) * eps(i, l);
      }
      eps(i, j) += eps_1(i, j) - tmp;
      tmp_eps_1(i, j) = eps(i, j);
      tmp = 0;
    }
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      for (int j1 = 0; j1 < pdim; j1++) {
        tmp_eps_2(i, j) += eps_2(i, j1) * b_D_dot(j1, j);
        eps(i, j) += eps_2(i, j1) * b_D_dot(j1, j);
      }
    }
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      eps_beta[i] += eps(i, j) * risk_cumsum[j];
    }
    eps_beta[i] /= tot_risk;
  }
  
  
  for (int j = 0; j < k; j++) {
    for (int l = 0; l < k; l++) {
      for (int i = 0; i < n; i++) {
        b_var_est(j, l) += eps(i, j) * eps(i, l);
      }
    }
  }
  
  
  
  for (int j = 0; j < k; j++) {
    for (int l = 0; l < j; l++) {
      B_D_se[j] += b_var_est(l, j);
    }
    for (int l = 0; l < j; l++) {
      B_D_se[j] += b_var_est(j, l);
    }
    B_D_se[j] += b_var_est(j, j);
    if (j > 0) {
      B_D_se[j] += B_D_se[j - 1];
    }
  }
  
  for (int j = 0; j < k; j++) {
    B_D_se[j] = sqrt(B_D_se[j]);
  }
  
  for (int i = 0; i < n; i++) {
    beta_se += eps_beta[i] * eps_beta[i];
  }
  
  beta_se = sqrt(beta_se);
  
  List by_prod = List::create(
    Named("Zc") = Zc,
    Named("dN") = dN,
    Named("risk_t") = risk_t,
    Named("indik_v") = indik_v,
    Named("tmp_eps_1") = tmp_eps_1,
    Named("tmp_eps_2") = tmp_eps_2,
    Named("b_numer_D") = b_numer_D,
    Named("b_dinom_D") = b_dinom_D,
    Named("B_intD") = B_intD,
    Named("eps_1") = eps_1,
    Named("eps_2") = eps_2,
    Named("eps") = eps,
    Named("H_dot") = H_dot,
    Named("b_D_dot") = b_D_dot,
    Named("b_var_est") = b_var_est,
    Named("eps_beta") = eps_beta);
  
  return List::create(
    Named("stime") = stime,
    Named("dB_D") = dB_D,
    Named("B_D") = B_D,
    Named("beta") = beta,
    Named("B_D_se") = B_D_se,
    Named("beta_se") = beta_se,
    Named("by_prod") = by_prod);
}
