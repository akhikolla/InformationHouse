#include <Rcpp.h>
using namespace Rcpp;
///' Gibbs sampling of the quantitative latent variables
///'
///' The function \code{CompRateYpartial} simulates from the Exponential full-conditionals of the quantitative latent variables for the Gibbs sampling of a Bayesian mixture of Plackett-Luce models.
///'
///' @param rate Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of rate parameters.
///' @param n_rank Numeric vector of length \eqn{N} with the number of items ranked by each sample unit.
///' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of posterior samples of the quantitative latent variables.
// [[Rcpp::export]]
  NumericMatrix SimYpsilon(NumericMatrix rate, NumericVector n_rank) {

    int N = rate.nrow();
    int K = rate.ncol();

    NumericMatrix out(N,K) ;

   int    s ;
   int    slot ;

        for( s=0 ; s<N ; s++){
           for( slot=0 ; slot<n_rank[s] ; slot++){
      out(s,slot) = R::rexp(1/rate(s,slot));
}
}

return out ; 

}
