#include <Rcpp.h>
using namespace Rcpp;

///' Binary matrix detailing the items ranked by each sample unit
///'
///' The function \code{umat} returns a binary matrix whose elements indicate which items has been ranked by each sample unit.
///'
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///'
///' @return Binary \eqn{N}\eqn{\times}{x}\eqn{K} matrix indicating whether sample unit \eqn{s} ranked item \eqn{i}.
// [[Rcpp::export]]
  NumericMatrix umat(NumericMatrix pi_inv) {

    int N = pi_inv.nrow();
    int K = pi_inv.ncol();

    NumericMatrix out(N, K) ;

   int    s ;
   int    slot ;

   for(s=0; s<N; s++){
     for(slot=0; slot<K; slot++){
     if(pi_inv(s,slot)>0){
     out(s,(int) pi_inv(s,slot)-1) = 1.0 ;
}
}
}

   return out;
  }
