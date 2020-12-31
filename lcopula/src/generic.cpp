// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
using namespace Rcpp;

//' Rcpp function to compute sum of rows and couple according to alphavec
//'
//' @param mc a matrix of Archimedean vectors to share of dimension \code{sum(alphavec)}
//' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
//' 
//' @keywords internal
//' 
//' @return a matrix of dimension \code{n} by \code{length(alphavec)} with combined data
// [[Rcpp::export(.marginCombo)]]
NumericMatrix marginCombo(NumericVector alphavec, NumericMatrix mc){
  int dim = alphavec.size();
  NumericVector cum_alpha=NumericVector(dim+1);
  for(int l=1;l<(dim+1);l++){
    cum_alpha[l]=cum_alpha[l-1]+alphavec[l-1];
  }
  NumericMatrix combined = NumericMatrix(mc.nrow(),dim);
  for(int i=0;i<dim;i++){
    for(int n=0;n<mc.nrow();n++){
      combined(n,i)=0;
      for(int j=cum_alpha[i];j<cum_alpha[i+1];j++){
        combined(n,i)+=mc(n,j);
      }
    }
  }
  return(combined);
}
