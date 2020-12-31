#include <Rcpp.h>
using namespace Rcpp;
///' Rate parameter computation for the Gamma full-conditionals of the support parameters
///'
///' The function \code{CompRateP} computes the rate parameters of the Gamma full-conditionals of the support parameters for the Gibbs sampling of a Bayesian mixture of Plackett-Luce models.
///'
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @param Y Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of the quantitative latent variables.
///' @param z Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of binary component memberships.
///' @param u_bin Binary \eqn{N}\eqn{\times}{x}\eqn{K} matrix indicating whether sample unit \eqn{s} ranked item \eqn{i}.
///' @param n_rank Numeric vector of length \eqn{N} with the number of items ranked by each sample unit.
///' @param rate0 Numeric vector of \eqn{G} rate hyperparameters.
///' @return Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of rate parameters.
// [[Rcpp::export]]
  NumericMatrix CompRateP(NumericMatrix pi_inv, NumericMatrix Y, NumericMatrix z, NumericMatrix u_bin, IntegerVector n_rank, NumericVector rate0) {

    int N = pi_inv.nrow() ;
    int K = pi_inv.ncol() ;
    int G = z.ncol() ;

    NumericMatrix den(G,K) ;    
    NumericMatrix out(G,K) ;    
    NumericMatrix pdenom(N,K) ;    
    NumericVector temp(N) ;    

    int    s ;
    int    slot ;
    int    item ;
    int    group ;
    double    availablenext ;

    for( group=0 ; group<G ; group++){

      for( item=0 ; item<K ; item++){

// compute colsums of z_hat multiplied by u matrix -> gamma.hat

       den(group,item) = 0.0 ;
       
         for( s=0 ; s<N ; s++){
        
           temp[s] = 0.0 ; 

          availablenext = 1.0 ;

            for( slot=0 ; slot<n_rank[s] ; slot++){
              if(availablenext == 1.0){        
                   temp[s] = temp[s] + Y(s,slot) ; 
                 }
                    if( pi_inv(s,slot) == ((double) (item+1)) ){ 
                       availablenext = 0.0 ;
                    }
              }            

       den(group,item) = den(group,item) + z(s,group)*temp[s] ;
}      
       out(group,item) = (den(group,item)+rate0[group]) ;
}
}

return out ; 

}
