# include <Rcpp.h>
# include "fittingmeasure.h"
using namespace Rcpp;
///' Conditional Chi-squared index relying on paired comparisons for observed data
///'
///'
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
///' @param weights Numeric vector of the \eqn{G} mixture weights.
///' @return Conditional Chi-squared index value.
// [[Rcpp::export]]
double chisqmeasureobscond(Rcpp::IntegerMatrix pi_inv, Rcpp::NumericMatrix p, Rcpp::NumericVector weights){

    double f = 0.0 ;
       int N = pi_inv.nrow();
       int K = pi_inv.ncol();
       int G = p.nrow();
       int slot ;
       int slot2 ;
       int g ;
       int s ;
       int hmr ;
       int temphmr ;

       Rcpp::IntegerMatrix tau_mat(K,K);
       Rcpp::NumericMatrix tau_star_mat(K,K);
       Rcpp::IntegerMatrix T_mat(K,K);
       Rcpp::IntegerMatrix temp_pi_inv(N,K);

       Rcpp::NumericVector marg_p(K); 

       for(slot=0; slot<K; slot++){
          for(g=0; g<G; g++){
	 marg_p[slot] = marg_p[slot] + weights[g]*p(g,slot) ;
          }
       }


       for(hmr=0; hmr<K; hmr++){      

         
    for(s=0; s<N; s++){
      temphmr=0 ;
      /* count how many ranked in each subject s */
       for(slot=0; slot<K; slot++){
         temp_pi_inv(s,slot)=0;        
	 if(pi_inv(s,slot)>0){
	   temphmr=temphmr+1 ;
           temp_pi_inv(s,slot)=pi_inv(s,slot);
	 }
       }
      /* if the items ranked by subject s are not hmr put all 0 */
	 if(temphmr!=(hmr+1)){         
	   for(slot=0; slot<K; slot++){
	     temp_pi_inv(s,slot)=0;
	   }
	 }
    }

       tau_mat = tau(temp_pi_inv); 
       for(slot=0; slot<K; slot++){
         for(slot2=0; slot2<slot; slot2++){
	   T_mat(slot,slot2) = tau_mat(slot,slot2) + tau_mat(slot2,slot);
	   T_mat(slot2,slot) = T_mat(slot,slot2) ;
	     tau_star_mat(slot,slot2) = ((double) T_mat(slot,slot2))*marg_p[slot]/(marg_p[slot]+marg_p[slot2]);
	   if(tau_star_mat(slot,slot2)>0.0){    
	   f = f + (((double) tau_mat(slot,slot2))-tau_star_mat(slot,slot2))*(((double) tau_mat(slot,slot2))-tau_star_mat(slot,slot2))/(tau_star_mat(slot,slot2));
	   	 }                               

	 }
       }
       }

       return f;  

}
