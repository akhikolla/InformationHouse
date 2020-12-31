# include <Rcpp.h>
# include "fittingmeasure.h"
using namespace Rcpp;
///' Computation top1 frequencies conditionally on the number of ranked items
///'
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///' @return Numeric \eqn{K}\eqn{\times}{x}\eqn{K} matrix of top1 frequencies.
// [[Rcpp::export]]

Rcpp::IntegerMatrix top1freq1dim(Rcpp::IntegerMatrix pi_inv) {
    int N = pi_inv.nrow();
    int K = pi_inv.ncol();

    Rcpp::IntegerMatrix freq(K,K);
    Rcpp::IntegerVector temp(K);

    int    s ;
    int    slot ;
    int    nranked ;

    for(s=0; s<N; s++){
      /* determine how many items have been ranked */

      nranked = 0;
      slot = 0 ;

      while(pi_inv(s,slot)>0 && slot<(K-1)){
	nranked = nranked + 1;
        slot = slot +1;
      }
      if(slot==(K-1) && pi_inv(s,slot)>0){
	nranked = nranked + 1;
      }

      /* compute/increment how many times item k is ranked as top1 */
      /* among the subjects who have ranked top-nranked items */

      freq(nranked-1,pi_inv(s,0)-1) = freq(nranked-1,pi_inv(s,0)-1) + 1;

    }

   return freq;
}


