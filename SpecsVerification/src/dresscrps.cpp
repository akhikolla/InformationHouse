#include <Rcpp.h>
using namespace Rcpp;

double _pnorm(double x) {
  return R::pnorm(x, 0.0, 1.0, 1, 0);
}
double _dnorm(double x) {
  return R::dnorm(x, 0.0, 1.0, 0);
}

//' Dress CRPS
//'
//' @param m vector of kernel means
//' @param s vector of kernel standard deviations
//' @param y observation 
//' @return crps
//' @export
// [[Rcpp::export]]
double dresscrps_cpp(NumericVector m, NumericVector s, double y) {
  int K = m.size();

  double sum1 = 0.0;
  
  for (int i = 0; i < K; ++i) {
    double zi = (m[i] - y) / s[i];
    sum1 += (m[i] - y) * (2.0 * _pnorm(zi) - 1) + 2.0 * s[i] * _dnorm(zi) - s[i] / K * 0.5 * M_2_SQRTPI;
  }
  sum1 /= K;

  double sum2 = 0.0;
  for (int i = 1; i < K; ++i) {
    for (int j = 0; j < i; ++j) {
      double eiej = m[j] - m[i];
      double sisj = sqrt(s[i]*s[i] + s[j]*s[j]);
      double argg = eiej / sisj;
      sum2 += eiej * (2.0 * _pnorm(argg) - 1) + 2.0 * sisj * _dnorm(argg); 
    }
  }
  sum2 = sum2 / K / K;

  return(sum1 - sum2);
}

