#include <Rcpp.h>
using namespace Rcpp;
///' M-step for the support parameters of a Bayesian mixture of Plackett-Luce models
///'
///' The function \code{UpPhetpartial} updates the support parameter estimates in the EM algorithm for MAP estimation of a Bayesian mixture of Plackett-Luce models.
///'
///' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
///' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders.
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @param u_bin Binary \eqn{N}\eqn{\times}{x}\eqn{K} matrix indicating whether sample unit \eqn{s} ranked item \eqn{i}.
///' @param z_hat Numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of posterior component membership probabilities.
///' @param shape0 Numeric \eqn{G}\eqn{\times}{x}{K} matrix of shape hyperparameters.
///' @param rate0 Numeric vector of \eqn{G} rate hyperparameters.
///' @param n_rank Numeric vector of length \eqn{N} with the number of items ranked by each sample unit.
///' @return Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of estimated component-specific support parameters.
// [[Rcpp::export]]
  NumericMatrix UpPhetpartial(NumericMatrix p, NumericMatrix ref_order, NumericMatrix pi_inv, NumericMatrix u_bin, NumericMatrix z_hat, NumericMatrix shape0, NumericVector rate0, IntegerVector n_rank) {

    int N = pi_inv.nrow() ;
    int K = pi_inv.ncol() ;
    int G = p.nrow() ;

    NumericMatrix num(G,K) ;    
    NumericMatrix den(G,K) ;    
    NumericMatrix out(G,K) ;    
    NumericMatrix pdenom(N,K) ;    
    NumericVector temp(N) ;    

    int    s ;
    int    slot ;
    int    slot2 ;
    int    item ;
    int    group ;
    double    availablenext ;
    double    temp_pdenom ;

    for( group=0 ; group<G ; group++){
      for( item=0 ; item<K ; item++){

// compute colsums of z_hat multiplied by u matrix -> gamma.hat

       den(group,item) = 0.0 ;
       num(group,item) = 0.0 ;

        for(s=0 ; s<N ; s++){
               
       num(group,item) = num(group,item) + z_hat(s,group)*u_bin(s,item) ;

       temp[s] = 0.0 ;
       temp_pdenom = 0.0;

/* compute the INITIAL denominator which is always needed for the first slot */
/* then update the denominator removing the support parameter of the item  */
/* selected at the current slot  */

       for( slot2=0 ; slot2<K ; slot2++){
           temp_pdenom = temp_pdenom + p(group,slot2) ; 
          }

   availablenext = 1.0 ;

/* any item is always available next at the first step/slot  */

         for( slot=0 ; slot<n_rank[s] ; slot++){

if(availablenext == 1.0){

            temp[s] = temp[s] + 1/temp_pdenom ;

/* but when the item is selected at the current slot */
/* it is no longer available next  */

         if( pi_inv(s,slot) == ((double) (item+1)) ){ 
           availablenext = 0.0 ;
       }
            temp_pdenom = temp_pdenom - p(group,((int) pi_inv(s,slot))-1) ;

         }
}

       den(group,item) = den(group,item) + z_hat(s,group)*temp[s] ;

}

       if(den(group,item)<=0){
       den(group,item) = 0.000000000001 ;
}

       out(group,item) = (num(group,item)+shape0(group,item)-1.0)/(den(group,item)+rate0[group]) ;

       if(out(group,item)<=0){
       out(group,item) = 0.000000000001 ;
}
}
}

return out ; 

}

