#include <Rcpp.h>
using namespace Rcpp;

//' CRPS for ensemble forecasts (C++ implementation)
//'
//' @param ens Ensemble members as columns of a matrix
//' @param obs The verifying observations
//' @param R_new Size for ensemble adjustment
//'
//' @return vector of crps values
//' @export
// [[Rcpp::export]]

NumericVector enscrps_cpp(NumericMatrix ens, NumericVector obs, double R_new) {

  int n = obs.size();
  int R_max = ens.ncol();

  NumericVector R(n);
  NumericVector sad_ens_half(n);
  NumericVector sad_obs(n);

  NumericVector crps(n);

  // calculate half sum of absolute difference between ensemble members:
  //   sad_ens_half = apply(ens, 1, function(x) sum(dist(x), na.rm=TRUE))
  // and sum of absolute differences between members and obs:
  //   sad_obs = rowSums(abs(ens - obs), na.rm=TRUE)
  // and vector of ensemble:
  //   R = rowSums(is.finite(ens))
  for (int i = 0; i < n; i++) {

    double the_obs = obs(i);

    NumericVector the_ens = ens(i,_);
    the_ens.sort(); // NAs are last in order

    sad_obs(i) = 0;
    double sum_xj = 0;
    double sum_jxj = 0;

    int j = 0;
    while (j < R_max) {
      if (R_IsNA(the_ens[j])) { 
        break;
      } 
      sad_obs(i) += fabs(the_ens(j) - the_obs);
      sum_xj     += the_ens(j);
      sum_jxj    += (j+1) * the_ens(j);
      R(i)       += 1;
      j          += 1;

    }
    sad_ens_half[i] = 2.0 * sum_jxj - (R(i)+1) * sum_xj;
  }

  // calculate the ensemble-adjusted crps 

  if (R_IsNA(R_new)) { 

    // no adjustment, default case, we catch this first because we assume that
    // this is the most common usage of the function
    for (int i = 0; i < n; i++) {
      crps(i) = sad_obs(i) / R(i) - sad_ens_half(i) / (R(i) * R(i));
    }

  } else if (R_new > 1) { // user wants ensemble-adjustment 

    // calculate the ensemble adjusted crps; returns the fair crps if
    // R.new==Inf, returns the unadjusted crps as above for R.new==R; returns NA
    // whenever R==1 (which is desired)
    for (int i = 0; i < n; i++) {
      crps(i) = sad_obs(i) / R(i) - sad_ens_half(i) / (R(i)*(R(i)-1)) * (1 - 1/R_new);
    }

  } else if (R_new == 1) {

    // probably a rare case; essentially the same as bove but we account for
    // the fact that (1-1/R.new)=0 so the entire sad.ens.half term vanishes;
    // also: for R.new == R == 1 the above would return NA, because one-member
    // ensembles can't be adjusted to new ensemble sizes; but technically we
    // have R == R.new, so no ensemble-adjustment was required, so the
    // unadjusted crps is returned, which is sad.obs/R as well for R==1
    for (int i = 0; i < n; i++) {
      crps(i) = sad_obs(i) / R(i);
    }

  } else {

    // all remaining cases such as R.new <= 0
    for (int i = 0; i < n; i++) {
      crps(i) = NA_REAL;
    }

  }

  return crps;

}

