#include <Rcpp.h>
using namespace Rcpp;
///' Compute paired comparison matrix for a partial ordering dataset
///'
///' @param pi_inv Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings.
///'
///' @return The \eqn{K}\eqn{\times}{x}\eqn{K} paired comparison matrix: number of times that item \eqn{i} is preferred to item \eqn{i'}.
// [[Rcpp::export]]
Rcpp::IntegerMatrix tau(Rcpp::IntegerMatrix pi_inv) {
    int N = pi_inv.nrow();
    int K = pi_inv.ncol();
    
    Rcpp::IntegerMatrix tau(K,K);
    Rcpp::IntegerVector starttemp(K);
    Rcpp::IntegerVector temp(K);
    
    int    s ;
    int    slot ;
    int    slot2 ;
    int    k ;
    
    for(s=0; s<N; s++){
        
        for(k=0; k<K; k++){
            starttemp[k] = 1;
            temp[k] = starttemp[k];
        }
        
        for(slot=0; slot<(K-1); slot++){
            if(pi_inv(s,slot)>0){
                
                starttemp[pi_inv(s,slot)-1] = 0 ;
                
                for(k=0; k<K; k++){
                    temp[k]=starttemp[k] ;
                }
                
                slot2 = slot+1 ;
                
                while(slot2<K){
                    if(pi_inv(s,slot)>0 && pi_inv(s,slot2)>0){
                        tau(pi_inv(s,slot)-1, pi_inv(s,slot2)-1) = tau(pi_inv(s,slot)-1, pi_inv(s,slot2)-1)+1 ;
                        temp[pi_inv(s,slot2)-1]=0;
                    }
                    slot2 = slot2+1 ;
                }
                if(pi_inv(s,slot)>0){
                    for(k=0; k<K; k++){
                        tau(pi_inv(s,slot)-1,k) = tau(pi_inv(s,slot)-1,k) + temp[k] ;
                    }
                }
            }
        }
    }
    return tau;
}

