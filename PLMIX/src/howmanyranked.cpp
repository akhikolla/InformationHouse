#include <Rcpp.h>
using namespace Rcpp;

///' Counts how many items are ranked in a partial ranking matrix
///' 
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @return Numeric vector of length \eqn{N} with the number of items ranked by each sample unit
// [[Rcpp::export]]
  IntegerVector howmanyranked(NumericMatrix pi_inv) {
    int N = pi_inv.nrow();
    int K = pi_inv.ncol();
    IntegerVector out(N) ;

   int    s ;
   int    slot;

   for(s=0; s<N; s++){
   for(slot=0; slot<K; slot++){
   if(pi_inv(s,slot)>0){
   out[s]=out[s]+1;
}
}
}
   return out;
}
