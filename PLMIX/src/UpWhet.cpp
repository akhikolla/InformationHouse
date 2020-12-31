#include <Rcpp.h>
using namespace Rcpp;
///' M-step for the weights of a Bayesian mixture of Plackett-Luce models
///'
///' The function \code{UpWhet} updates the mixture weight estimates in the EM algorithm for MAP estimation of a Bayesian mixture of Plackett-Luce models
///'
///' @param z_hat Numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of posterior component membership probabilities.
///' @param alpha0 Numeric vector of \eqn{G} Dirichlet hyperparameters.
///' @return Numeric vector of the \eqn{G} estimated mixture weights.
// [[Rcpp::export]]
  NumericVector UpWhet(NumericMatrix z_hat, NumericVector alpha0) {

    int N = z_hat.nrow();
    int G = z_hat.ncol();
    NumericVector out(G);

    int s ;
    int group ; 
    double sum=0.0 ;

    for( group=0 ; group<G ; group++){
      for( s=0 ; s<N ; s++){
       out[group] = out[group] + z_hat(s,group) ;
      }
if(ISNAN(out[group])){
out[group] = 0.0000000000000000000001;
}
    sum = sum + out[group] ;
}

    for( group=0 ; group<G ; group++){
       out[group] = out[group]/sum;
}

    return out ;
}
