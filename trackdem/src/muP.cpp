// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector muP(NumericVector m, NumericVector id,
                  NumericVector cm1, NumericVector cm2, NumericVector cm3, 
                  IntegerVector d) {
  using namespace Rcpp;
  using namespace arma;
 
  NumericVector dim(d);
  NumericVector i(id);
  NumericMatrix mat(m);
  NumericMatrix cmat1(cm1);
  NumericMatrix cmat2(cm2);
  NumericMatrix cmat3(cm3);

  arma::mat y1(mat.begin(),dim[0],dim[1],false);
 
  NumericMatrix mu(i.size(),3);
  
  for (int j = 0; j < i.size(); j++) { 
    // find coordinates
    arma::uvec ind = arma::find(y1 == i[j]);
  
      // calculate mean values
      for (int k = 0; k < ind.size(); k++) {
        mu(j,0) = mu(j,0) + cmat1[ind[k]];
        mu(j,1) = mu(j,1) + cmat2[ind[k]];
        mu(j,2) = mu(j,2) + cmat3[ind[k]];
      }
    // mean pixel values
    mu(j,0) = mu(j,0) / ind.size();
    mu(j,1) = mu(j,1) / ind.size();
    mu(j,2) = mu(j,2) / ind.size();
  
}
  return(wrap(mu));
}
