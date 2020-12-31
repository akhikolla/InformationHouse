#include <Rcpp.h>
using namespace Rcpp;

template <class T>
inline double do_yFNVOYNjnZFHTeDvjESk( T& x ) {
  int n = x.size();
  LogicalVector out(n);
  double k=0;
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
    if ( out[i] == FALSE ){
      ++k;
    }
  }
  return k;
}

NumericVector row_yFNVOYNjnZFHTeDvjESks( NumericMatrix& X ) {

  int nRows = X.nrow();
  NumericVector out = no_init(nRows);

  for( int i=0; i < nRows; i++ ) {
    NumericMatrix::Row tmp = X(i, _);
    out[i] = do_yFNVOYNjnZFHTeDvjESk( tmp );
  }

  return out;

}

NumericVector col_yFNVOYNjnZFHTeDvjESks( NumericMatrix& X ) {

  int nCols = X.ncol();
  NumericVector out = no_init(nCols);

  for( int j=0; j < nCols; j++ ) {
    NumericMatrix::Column tmp = X(_, j);
    out[j] = do_yFNVOYNjnZFHTeDvjESk( tmp );
  }

  return out;

}

// [[Rcpp::export]]
NumericVector fast_apply_nb_not_na( NumericMatrix X, int dim ) {

  NumericVector a;
  a=NULL;
  if( dim == 1 ) {
    a=row_yFNVOYNjnZFHTeDvjESks(X);
  } else if( dim == 2 ) {
    a=col_yFNVOYNjnZFHTeDvjESks(X);
  }

  return a;
}

