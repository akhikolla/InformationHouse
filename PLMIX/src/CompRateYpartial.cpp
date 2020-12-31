#include <Rcpp.h>
using namespace Rcpp;
///' Rate parameter computation for the Exponential full-conditionals of the quantitative latent variables
///'
///' The function \code{CompRateYpartial} computes the rate parameters of the Exponential full-conditionals of the quantitative latent variables for the Gibbs sampling of a Bayesian mixture of Plackett-Luce models.
///'
///' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders.
///' @param z Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of binary component memberships.
///' @param n_rank Numeric vector of length \eqn{N} with the number of items ranked by each sample unit.
///' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of rate parameters.
// [[Rcpp::export]]
  NumericMatrix CompRateYpartial(NumericMatrix p, NumericMatrix pi_inv, NumericVector ref_order, NumericMatrix z, NumericVector n_rank) {

    int N = pi_inv.nrow();
    int K = pi_inv.ncol();
    int G = p.nrow();

    NumericMatrix out(N,K) ;
    NumericMatrix pdenom(N,K) ;

   int    s ;
   int    group ;
   int    slot ;
   int    slot3 ;
   NumericVector initialsum(G) ;
   double tempsum ;


        for( group=0 ; group<G ; group++){
           for( slot3=0 ; slot3<K ; slot3++){
              initialsum[group] = initialsum[group] + p(group,slot3) ; 
            }
        }

  for( s=0 ; s<N ; s++){
      for( group=0 ; group<G ; group++){
         if(z(s,group)>0){
           tempsum = initialsum[group] ;
           for( slot=0 ; slot<n_rank[s] ; slot++){            
                  out(s,slot) = tempsum ;
                  tempsum = tempsum - p(group,((int) pi_inv(s,slot))-1) ; 
            }
         }
      }
  }

return out ; 

}
