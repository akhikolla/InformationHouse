# include <Rcpp.h>
# include "fittingmeasure.h"
using namespace Rcpp;
///' Chi-squared index relying on paired comparisons for observed data
///'
///'
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
///' @param weights Numeric vector of the \eqn{G} mixture weights.
///' @return Chi-squared index value.
// [[Rcpp::export]]
double chisqmeasureobs(Rcpp::IntegerMatrix pi_inv, Rcpp::NumericMatrix p, Rcpp::NumericVector weights){

    double f = 0.0 ;
       int K = pi_inv.ncol();
       int G = p.nrow();
       int slot ;
       int slot2 ;
       int g ;

       Rcpp::IntegerMatrix tau_mat(K,K);
       Rcpp::NumericMatrix tau_star_mat(K,K);
       Rcpp::IntegerMatrix T_mat(K,K);

       Rcpp::NumericVector marg_p(K); 


       for(slot=0; slot<K; slot++){
          for(g=0; g<G; g++){
	 marg_p[slot] = marg_p[slot] + weights[g]*p(g,slot) ;
          }
       }


       tau_mat = tau(pi_inv); 
       for(slot=0; slot<K; slot++){
         for(slot2=0; slot2<slot; slot2++){
	   T_mat(slot,slot2) = tau_mat(slot,slot2) + tau_mat(slot2,slot);
	   T_mat(slot2,slot) = T_mat(slot,slot2) ;
	   tau_star_mat(slot,slot2) = ((double) T_mat(slot,slot2))*marg_p[slot]/(marg_p[slot]+marg_p[slot2]);
	   f = f + (((double) tau_mat(slot,slot2))-tau_star_mat(slot,slot2))*(((double) tau_mat(slot,slot2))-tau_star_mat(slot,slot2))/(tau_star_mat(slot,slot2));
	 }
       }



       return f;  

}

