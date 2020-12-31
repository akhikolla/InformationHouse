#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector find_corres(NumericMatrix X0, NumericMatrix X) {

  // @param X0 matrix of unique designs
  // @param X matrix of all designs
  // @return vector associating rows of X with those of X0
  
  
  IntegerVector corres(X.nrow());
  
  bool tmp;
  int nX0r = X0.nrow();
  int nX0c = X0.ncol();
  
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < nX0r; j++){
      tmp = true;
      for(int k = 0; k < nX0c; k++){
        if(X(i,k) != X0(j,k)){
          tmp = false;
          break;
        }
      }
      if(tmp){
        corres(i) = j + 1;
        break;
      }
    }
  }
  return(corres);
}


