#include <Rcpp.h>
using namespace Rcpp ;

template <int RTYPE>
Vector<RTYPE> na_filldown_template(const Vector<RTYPE>& x) {
  int n = x.size() ;
  int n_out = n ;
  int lastgood = 0;

  Vector<RTYPE> out(n_out) ;
  for( int i=0, j=0; i<n; i++){
    if( Vector<RTYPE>::is_na( x[i] ) )
      out[j++] = x[lastgood];
    else {
      out[j++] = x[i];
      lastgood=i;
    }
  }
  return out ;
}



//' Fill down NA with the last observed observation
//'
//' @description Fill down missing values with the latest non-missing value
//' @param x A vector
//' @return A vector or list with the NA's replaced by the last observed value.
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' a <- c(1:5, "Howdy", NA, NA, 2:3, NA)
//' filldown(a)
//' filldown(c(NA, NA, NA, 3:5))
//'
//' @export
// [[Rcpp::export]]
SEXP filldown( SEXP x ){
  RCPP_RETURN_VECTOR( na_filldown_template, x ) ;
}



