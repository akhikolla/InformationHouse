#include "mosum_util.h"

//' equivalent to rollsum(x, k=G, fill=NA, align="left") in the package zoo, 
//' but optimized for speed
//' @keywords internal
// [[Rcpp::export]]
NumericVector rolling_sum(const NumericVector &x, unsigned G) {
  const unsigned n = x.length();
  NumericVector res(n, NA_REAL);
  double currentSum = 0.0;
  for (unsigned j=0; j<G; ++j) {
    currentSum += x[j];
  }
  res[0] = currentSum;
  for (unsigned j=1; j<n-G+1; ++j) {
    currentSum += x[j+G-1];
    currentSum -= x[j-1];
    res[j] = currentSum;
  }
  return res;
}

//' extract changepoints from candidates with eta criterion
//' @keywords internal
// [[Rcpp::export]]
IntegerVector eta_criterion_help(const IntegerVector &candidates, 
                                 const NumericVector &m_values,
                                 double eta, double G_left, double G_right) {
  const int n = m_values.length();
  IntegerVector res(0);
  const int left_length = std::floor(eta*G_left);
  const int right_length = std::floor(eta*G_right);
  for (int j=0; j<candidates.length(); ++j) {
    const int k_star = candidates[j];
    // Careful: k_star 1-indexed, vectors here are 0-indexed
    const double m_star = m_values[k_star-1]; 
    const int left_thresh = std::max(1, k_star-left_length) - 1; // see above
    const int right_thresh = std::min(n, k_star+right_length) - 1; // ""
    bool largest = true;
    for (int l=left_thresh; (l<=right_thresh) && largest; ++l) {
      if (m_values[l] > m_star) {
        largest = false;
      }
    }
    if (largest) {
      // k_star is maximum in its eta*G environment --> accept as changepoint
      res.push_back(k_star);
    }
  }
  return res;
}
