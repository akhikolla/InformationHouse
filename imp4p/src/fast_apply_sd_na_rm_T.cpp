#include <Rcpp.h>
using namespace Rcpp;

template <class T>
  inline double do_AXIlKpnLusRGZOyHYzIU( T& x ) {
    int n = x.size();
    LogicalVector out(n);
    double k=0;
    double l=0;
    for (int i = 0; i < n; ++i) {
      out[i] = NumericVector::is_na(x[i]);
      if ( out[i] == FALSE ){
        k += x[i];
        l++;
      }
    }
    double m=k/l;
    k=0;
    for (int i = 0; i < n; ++i) {
      out[i] = NumericVector::is_na(x[i]);
      if ( out[i] == FALSE ){
        k += (x[i]-m)*(x[i]-m);
      }
    }
    double res;
    if (l>1){res=sqrt(k/(l-1));}
    else{res=NA_REAL;}
    return res;
  }

NumericVector row_AXIlKpnLusRGZOyHYzIUs( NumericMatrix& X ) {

  int nRows = X.nrow();
  NumericVector out = no_init(nRows);

  for( int i=0; i < nRows; i++ ) {
    NumericMatrix::Row tmp = X(i, _);
    out[i] = do_AXIlKpnLusRGZOyHYzIU( tmp );
  }

  return out;

}

NumericVector col_AXIlKpnLusRGZOyHYzIUs( NumericMatrix& X ) {

  int nCols = X.ncol();
  NumericVector out = no_init(nCols);

  for( int j=0; j < nCols; j++ ) {
    NumericMatrix::Column tmp = X(_, j);
    out[j] = do_AXIlKpnLusRGZOyHYzIU( tmp );
  }

  return out;

}

// [[Rcpp::export]]
NumericVector fast_apply_sd_na_rm_T( NumericMatrix X, int dim ) {

  NumericVector a;
  a=NULL;
  if( dim == 1 ) {
    a=row_AXIlKpnLusRGZOyHYzIUs(X);
  } else if( dim == 2 ) {
    a=col_AXIlKpnLusRGZOyHYzIUs(X);
  }

  return a;
}

