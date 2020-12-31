# include <Rcpp.h>
# include "fittingmeasure.h"
using namespace Rcpp;
///' Random generation from a finite mixture of Plackett-Luce models and subsequent censoring
///'
///' Random generation from a finite mixture of Plackett-Luce models and subsequent censoring according to a given partial ordering matrix.
///'
///' @param N Number of sample units.
///' @param K Number of possible items.
///' @param G Number of mixture components.
///' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
///' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders.
///' @param weights Numeric vector of the \eqn{G} mixture weights.
///' @param rankingFormat Logical: whether the final simulated data should be expressed in ranking format.
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings (to replicate the observed missingness patterns).
///' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of simulated data (default is in ordering format) with the same missingness patterns of \code{pi_inv}.
// [[Rcpp::export]]
Rcpp::IntegerMatrix PLMIXsim(int N, int K, int G, Rcpp::NumericMatrix p, Rcpp::NumericMatrix ref_order,Rcpp::NumericVector weights, bool rankingFormat,  Rcpp::IntegerMatrix pi_inv){

Rcpp::IntegerMatrix out(N,K) ;
Rcpp::NumericVector rs(G) ;
Rcpp::NumericMatrix p_stand(G,K) ;
Rcpp::IntegerVector vwc(1) ;
Rcpp::IntegerVector temp(K) ;

      int g ;
      int s ;
      int slot ;
      int whichcomponent ;

     for(g=0; g<G; g++){
       for(slot=0; slot<K; slot++){
       rs[g] = rs[g] + p(g,slot) ;
   }
       for(slot=0; slot<K; slot++){
       p_stand(g,slot) =  p(g,slot)/rs[g];
       }
   }

        for(s=0; s<N; s++){
          vwc = quickintsample(G,1,weights);
          whichcomponent = vwc[0] ;
	  temp = quickintsample(K,K,p.row(whichcomponent-1));
       for(slot=0; slot<K; slot++){
	 out(s,ref_order(whichcomponent-1,slot)-1) = temp[slot] ;
	 if(pi_inv(s,ref_order(whichcomponent-1,slot)-1)==0){
	   out(s,ref_order(whichcomponent-1,slot)-1) = 0 ;
	 }
       }
	}

   return out;
}
