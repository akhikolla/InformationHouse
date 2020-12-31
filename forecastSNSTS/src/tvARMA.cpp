
#include <Rcpp.h>
using namespace Rcpp;

//################################################################################
//' Workhorse function for tvARMA time series generation
//'
//' More explanation!
//'
//' @name tvARMAcpp
//'
//' @param z a ...
//' @param x_int a ...
//' @param A ...
//' @param B a ...
//' @param Sigma a ...
//'
//' @return Returns a ...
//'
//################################################################################

// [[Rcpp::export]]
NumericVector tvARMAcpp( NumericVector z, NumericVector x_init, NumericMatrix A,
                      NumericMatrix B, NumericVector Sigma) {

   int T = A.nrow()-1; 
   NumericVector x( T );
   int p = A.ncol()-1;
   int q = B.ncol()-1;

   
   for (int i = 0; i < std::max(p,q); i++) {
      x(i) = x_init(i);
   }

   // Rcout << "The value is " << p << " " << q << std::endl;

   for (int i = std::max(p,q); i < T; i++)  {
    x(i) = Sigma(i+1) * z(q+i);
    if (q > 0) {
      for (int j = 1; j <= q; j++) {
        x(i) += B(i+1, j) * Sigma(i+1-j) * z(q+i-j);
      }
    }
    if (p > 0) {
      for (int j = 1; j <= p; j++) {
        x(i) += A(i+1, j) * x(i-j);
      }
    }
  }
   
  return x;

}
