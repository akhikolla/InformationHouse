#include <Rcpp.h>
using namespace Rcpp;
///' Multinomial probability computation for the Multinomial full-conditionals of the latent component memberships
///'
///' The function \code{CompProbZpartial} computes the multinomial probabilities of the Multinomial full-conditionals of the latent component memberships for the Gibbs sampling of a Bayesian mixture of Plackett-Luce models.
///'
///' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @param Y Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of the quantitative latent variables.
///' @param u_bin Binary \eqn{N}\eqn{\times}{x}\eqn{K} matrix indicating whether sample unit \eqn{s} ranked item \eqn{i}.
///' @param n_rank Numeric vector of length \eqn{N} with the number of items ranked by each sample unit.
///' @param omega Numeric vector of the \eqn{G} mixture weights.
///' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of multinomial probabilities.
// [[Rcpp::export]]
  NumericMatrix CompProbZpartial(NumericMatrix p, NumericMatrix pi_inv, NumericMatrix Y, NumericMatrix u_bin, IntegerVector n_rank, NumericVector omega) {

    int N = pi_inv.nrow() ;
    int K = pi_inv.ncol() ;
    int G = p.nrow() ;

    NumericMatrix out(N,G) ;    
    NumericVector temp(K) ;    

    int    s ;
    int    slot ;
    int    item ;
    int    group ;
    double    availablenext ;

    for( s=0 ; s<N ; s++ ){
      for( group=0 ; group<G ; group++ ){
        out(s,group) = 1.0 ;

        for( item=0 ; item<K ; item++){

          availablenext = 0.0 ;	
          temp[item] = 0.0 ;

          availablenext = 1.0 ;
            for( slot=0 ; slot<n_rank[s] ; slot++){
              if(availablenext == 1.0){        
                   temp[item] = temp[item] + Y(s,slot) ; 
                 }
                    if( pi_inv(s,slot) == ((double) (item+1)) ){ 
                       availablenext = 0.0 ;
                    }
              }

          if (u_bin(s,item)>0.0){
          out(s,group) = out(s,group)*p(group,item) ;
            }
            out(s,group) = out(s,group)*exp(-p(group,item)*temp[item]);
}
            out(s,group) = out(s,group) * omega[group];
}     
}

return out ; 

}
