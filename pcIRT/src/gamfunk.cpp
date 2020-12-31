#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gamfunk(NumericVector epsmat, NumericMatrix gammat)
{
            int len = epsmat.size();
            
            for(int i=1; i < len; ++i){
            for(int j=1; j < len; ++j){
            if(j >= i){
            if(i == j){
            gammat(i,j) = gammat(i-1, j-1)*epsmat(j);
            } else {
            gammat(i,j) = gammat(i, j-1) + gammat(i-1, j-1)*epsmat(j);
            }
            }
            }
            }
            return gammat(_,(len-1));
}

