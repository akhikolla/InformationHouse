#include <Rcpp.h>
#include <cmath>
#include "mosum_util.h"
using namespace Rcpp;

//' helping function for bootstrap (compute local means)
//' @keywords internal
// [[Rcpp::export]]
double mean_help(NumericVector x, int l, int r) {
  if (l>r) throw std::runtime_error("Expecting l<=r");
  double res = 0.0;
  for (int t=l; t<=r; ++t) {
    res += x[t];
  }
  res /= ((double)l - (double)r + 1);
  return res;
}

//' Compute bootstrapped mosum statistic and return maximum position thereof
//' @keywords internal
// [[Rcpp::export]]
int get_k_star(NumericVector x_star, // bootstrap replicated ts (of full length n)
           int k_hat, // estimated cpts pos
           int G_l, // detection bandwidth
           int G_r) {
  const int n = x_star.length();
  const int k_hat_ind = k_hat-1; // Recall that R is 1-indexed
  const int l = std::max(0, k_hat_ind-G_l+1);
  const int r = std::min(n-1, k_hat_ind+G_r);
  double max_val = -1.0;
  int max_pos = l-1;
  double current_val = -1.0;
  for (int t=l; t<=r; ++t) {
    // left boundary?
    if (t < G_l-1) {
      double t_val_help = 0.0;
      const double scaling = std::sqrt(((double)G_l + (double)G_r) / 
                                       ( (double)(t+1)*(double)(G_l+G_r-t-1) ));
      const double mean_l = mean_help(x_star, 0, G_l+G_r-1);
      for (int j=0; j<t; ++j) {
        t_val_help += (mean_l - x_star[j]);
      }
      current_val = std::abs(scaling * t_val_help);
    } // right boundary?
    else if (t >= n-G_r) {
      double t_val_help = 0.0;
      const double scaling = std::sqrt( ((double)G_l + (double)G_r) /
                                        ( (double)(n-t)*double(n-G_l-G_r+t) ));
      const double mean_r = mean_help(x_star, n-G_l-G_r, n-1);
      for (int j=n-G_l-G_r; j<n; ++j) {
        t_val_help += (mean_r - x_star[j]);
      }
      current_val = std::abs(scaling * t_val_help);
    } // middle?
    else {
      const double scaling = std::sqrt((double)G_l * (double)G_r / ((double)G_l + (double)G_r));
      const double mean_r = mean_help(x_star,t+1,t+G_r);
      const double mean_l = mean_help(x_star,t-G_l+1,t);
      current_val = std::abs(scaling * (mean_r-mean_l));
    }
    if (current_val > max_val) {
      max_val = current_val;
      max_pos = t;
    }
  }
  return max_pos+1; // Recall that R is 1-indexed
}

//' Obtain bootstrap replicate of time series
//' @keywords internal
// [[Rcpp::export]]
NumericVector bootstrapped_timeSeries(IntegerVector cpts,
                                      NumericVector x) {
  const int n = x.length();
  const int q = cpts.length();
  NumericVector x_star(n);
  for (int j=0; j<=q; ++j) {
    const int l = (j==0) ? 0 : (cpts[j-1]); // boundary for (stationary) residual part
    const int r = (j==q) ? n-1 : cpts[j]-1; // boundary for (stationary) residual part
    const int N = r-l+1;
    // Rcout << j << ", [" << l << "," << r << "]" << std::endl;
    const NumericVector x_local = x[Range(l, r)];
    x_star[Range(l, r)] = sample(x_local, N, true);
  }
  return x_star;
}

//' Helping function to get bootstrap replicates of change point estimates
//' @keywords internal
// [[Rcpp::export]]
List cpts_bootstrap_help(IntegerMatrix cpts_info,
                             NumericVector x,
                             int N_reps) {
  const IntegerVector cpts = cpts_info(_, 0);
  const int q = cpts.length();
  const int n = x.length();
  
  IntegerMatrix k_star(N_reps, q); // bootstrapped changepoints
  IntegerMatrix k_star1(N_reps, q); // bootstrapped differences (k^* - \hat k)
  NumericMatrix k_star2(N_reps, q); // bootstrapped rescaled differences d^2/sigma^2 (k^* - \hat k)
  
  NumericVector d_hat(q); // estimated jump heights
  NumericVector sigma2_hat(q); // estimated noise variance (pooled from left and right stationary part)
  for (int j=0; j<q; ++j) {
    const int m_pos = cpts[j]-1;
    const int l_pos = (j==0) ? 0 : (cpts[j-1]);
    const int r_pos = (j==q-1) ? n-1 : cpts[j+1]-1;
    const NumericVector x_l = x[Range(l_pos, m_pos)];
    const NumericVector x_r = x[Range(m_pos+1, r_pos)];
    d_hat[j] = mean(x_r) - mean(x_l);
    
    const double denominator = std::max(1.0, (double)(r_pos-l_pos-2));
    const double tau_l = var(x_l) * (double)(x_l.length()-1);
    const double tau_r = var(x_r) * (double)(x_r.length()-1);
    sigma2_hat[j] = (tau_l + tau_r) / denominator;
  }
  
  for (int iboot=0; iboot<N_reps; ++iboot) {
    // 
    // semi-local bootstrap (resample all but recompute locally)
    // note: redundant, may be sped up in future revisions
    //
    NumericVector x_star = bootstrapped_timeSeries(cpts, x);

    for (int j=0; j<q; ++j) {
      const int G_l = cpts_info(j, 1);
      const int G_r = cpts_info(j, 2);
      k_star(iboot,j) = get_k_star(x_star, cpts[j], G_l, G_r);
      k_star1(iboot,j) = k_star(iboot,j)-cpts[j];
      k_star2(iboot,j) = (double)(k_star(iboot,j)-cpts[j]) *
        std::pow(d_hat[j],2.0) / sigma2_hat[j];
    }
  }
  List res;
  res["k_star"] = k_star; // bootstrapped changepoints (raw)
  res["k_star1"] = k_star1; // bootstrapped differences (k_j^*-\hat k_j)
  res["k_star2"] = k_star2; // bootstrapped rescaled differences \hat d_j^2/ \hat \sigma_j^2(k_j^*-\hat k_j)
  res["d_hat"] = d_hat;
  res["sigma2_hat"] = sigma2_hat;
  return res;
}
