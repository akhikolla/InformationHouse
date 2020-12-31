#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate AUC difference `AUC(fcst,obs) - AUC(fcst_ref, obs)` of two forecasts for the same observations, and the sampling standard deviation of the AUC difference (Internal C++ implementation)
//' 
//' @param fcst numeric vector of forecasts (NAs are not allowed)
//' @param fcst_ref numeric vector of reference forecasts (NAs are not allowed)
//' @param obs vector of binary observations (obs[t] evaluates to TRUE if event happens at instance t, to FALSE otherwise) 
//' @return AUC values, their sampling standard deviations, the AUC difference, and their sampling standard deviations
//' @seealso Auc AucDiff
//' @export
// [[Rcpp::export]]
std::vector<double> aucdiff_cpp(std::vector<double> fcst, std::vector<double> fcst_ref, std::vector<bool> obs) {

  int L = obs.size();

  // calculate order vectors of fcst and fcst_ref
  arma::uvec i_ord = arma::sort_index(arma::vec(fcst), "ascend");
  arma::uvec i_ord_ref = arma::sort_index(arma::vec(fcst_ref), "ascend");

  int n,m,nn,mm,i,j,jp1,k;
  double sumV, sumV2, sumW, sumW2, auc, v, w, sd_auc,
         sumV_ref, sumV2_ref, sumW_ref, sumW2_ref, auc_ref, v_ref, w_ref, sd_auc_ref,
         sumV_V_ref, sumW_W_ref, v12, w12, auc_diff, sd_auc_diff;

  // vectors or length L to save V_i and W_j values 
  std::vector<double> VW(L, 0);
  std::vector<double> VW_ref(L, 0);

  // temporary vectors to save indices while looping over duplicates
  std::vector<double> V_inds(0);
  std::vector<double> W_inds(0);

  // calculate V and W for AUC(fcst, obs)
  n = m = i = 0;
  while (i < L) {
    nn = mm = 0;
    V_inds.clear();
    W_inds.clear();
    while (1) {
      j = i_ord[i];
      if (obs.at(j)) {
        mm++;
        V_inds.push_back(j);
      } else {
        nn++;
        W_inds.push_back(j);
      }
      if (i == L-1) {
        break;
      } 
      jp1 = i_ord[i+1];
      if (fcst.at(j) != fcst.at(jp1)) {
        break;
      } 
      i++;
    }
    for (k = 0; k < mm; k++) {
      VW.at(V_inds.at(k)) = n + nn/2.0;
    }
    for (k = 0; k < nn; k++) {
      VW.at(W_inds.at(k)) = m + mm/2.0;
    }
    n += nn;
    m += mm;
    i++;
  }

  // calculate V_ref and W_ref for AUC(fcst_ref, obs)
  n = m = i = 0;
  while (1) {
    nn = mm = 0;
    V_inds.clear();
    W_inds.clear();
    while (1) {
      j = i_ord_ref[i];
      if (obs.at(j)) {
        mm++;
        V_inds.push_back(j);
      } else {
        nn++;
        W_inds.push_back(j);
      }
      if (i == L-1) {
        break;
      } 
      jp1 = i_ord_ref[i+1];
      if (fcst_ref.at(j) != fcst_ref.at(jp1)) {
        break;
      } 
      i++;
    }
    for (k = 0; k < mm; k++) {
      VW_ref.at(V_inds.at(k)) = n + nn/2.0;
    }
    for (k = 0; k < nn; k++) {
      VW_ref.at(W_inds.at(k)) = m + mm/2.0;
    }
    n += nn;
    m += mm;
    i++;
    if (i >= L) {
      break;
    }
  }


  // calculate sums, sums-of-squares, and product sums to calculate variances
  // and covariances of V and W
  sumV = 0.0;
  sumV2 = 0.0;
  sumW = 0.0;
  sumW2 = 0.0;
  sumV_ref = 0.0;
  sumV2_ref = 0.0;
  sumW_ref = 0.0;
  sumW2_ref = 0.0;
  sumV_V_ref = 0.0;
  sumW_W_ref = 0.0;

  for (i = 0; i < L; i++) {
    if (obs.at(i)) {
      double V_ = VW.at(i) / n;
      double Vref_ = VW_ref.at(i) / n;
      sumV += V_;
      sumV2 += V_*V_;
      sumV_ref += Vref_;
      sumV2_ref += Vref_*Vref_;
      sumV_V_ref += V_ * Vref_;
    } else {
      double W_ = VW.at(i) / m;
      double Wref_ = VW_ref.at(i) / m;
      sumW += W_;
      sumW2 += W_*W_;
      sumW_ref += Wref_;
      sumW2_ref += Wref_*Wref_;
      sumW_W_ref += W_ * Wref_;
    }
  }

  // calculate quantities defined in Delong et al (1988)
  auc = sumV / m;
  auc_ref = sumV_ref /  m;
  v = (sumV2 - sumV * sumV / m) / (m-1);
  v_ref = (sumV2_ref - sumV_ref * sumV_ref / m) / (m-1);
  w = (sumW2 - sumW * sumW / n) / (n-1);
  w_ref = (sumW2_ref - sumW_ref * sumW_ref / n) / (n-1);
  v12 = (sumV_V_ref - sumV * sumV_ref / m) / (m-1);
  w12 = (sumW_W_ref - sumW * sumW_ref / n) / (n-1);
  
  
  // calculate AUC variances and covariances
  sd_auc = sqrt(v / m + w / n);
  sd_auc_ref = sqrt(v_ref / m + w_ref / n);
  auc_diff = auc - auc_ref;
  sd_auc_diff = sqrt((v + v_ref - 2 * v12) / m + (w + w_ref - 2 * w12) / n);


  // return
  double tmp[] = {auc, sd_auc, auc_ref, sd_auc_ref, auc_diff, sd_auc_diff};
  std::vector<double> ret(tmp, tmp + 6);
  return ret;
}


