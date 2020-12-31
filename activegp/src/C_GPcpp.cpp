# include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]
using namespace Rcpp;

NumericMatrix W_kappa_ij(NumericMatrix design, NumericVector theta, int i1, int i2, int ct);
double Ikk_cpp(double a, double b, double t, int ct);
double w_ii_cpp(double a, double b, double t, int ct);
double w_ij_cpp(double a, double b, double t, int ct);

//' Computes Int(kappa_i(X, design) . kappa_j(design, X)). This function is preferred for initialization
//' @title Covariance of kernel computations
//' @param design matrix of design points
//' @param theta lengthscales
//' @param Ki The inverse covariance matrix
//' @param Kir The inverse covariance matrix times the response.
//' @param ct Covariance type, 1 means Gaussian, 2 means Matern 3/2, 3 means Matern 5/2
//' @return The matrix representing the result of the integration.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix C_GP_cpp(NumericMatrix design, NumericVector response, NumericVector theta, NumericMatrix Ki, int ct){
  int n = design.nrow();
  int d = design.ncol();
  
  NumericVector Kir =  Ki * response;

  NumericMatrix C(d, d);
  
  double coef;
  double w;
  double wnew;
  // double wold;
  double ires;
  NumericVector wvec(d);
  
  // Compute main terms.
  for (int n1 = 0; n1 < n; n1++) {
    for (int n2 = 0; n2 < n; n2++) {
      coef = (Kir(n1) * Kir(n2) - Ki(n1,n2));
      
      w = 1.0;
      for (int d0 = 0; d0 < d; d0++) {
        ires = Ikk_cpp(design(n1, d0), design(n2, d0), theta(d0), ct);
        w *= ires;
        wvec(d0) = ires;
      }
      
      for (int d1 = 0; d1 < d; d1++) {
        for (int d2 = d1; d2 < d; d2++) {
          if (d1 == d2) {
            wnew = w_ii_cpp(design(n1, d1), design(n2, d1), theta(d1), ct);
            C(d1,d2) += coef * w * wnew / wvec[d1];
          } else {
            wnew = w_ij_cpp(design(n1, d1), design(n2, d1), theta(d1), ct) * w_ij_cpp(design(n2, d2), design(n1, d2), theta(d2), ct);
            C(d2, d1) = C(d1,d2) += coef * w * wnew / (wvec[d1] * wvec[d2]);
          }
        }
      }
    }
  }
  
  // Augment diag
  for (int d0 = 0; d0 < d; d0++) {
    C(d0,d0) += 1 / pow(theta(d0), 2);
  }
  
  return(C);
}

