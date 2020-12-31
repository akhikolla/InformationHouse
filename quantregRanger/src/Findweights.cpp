#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector Findweights(NumericVector ONv, NumericVector NNv, NumericVector WEv, int nobs, int nnew, int ntree, double thres, NumericVector counti, int normalise){
  
  int k,i,j;
  double absol;
  int meancount;
  
  for (k=1; k<=ntree; ++k){
    for (i=1; i<=nnew; ++i){
      meancount = 0;
      
      for (j=1; j<=nobs; ++j){
        absol = ONv[j+(k-1)* (nobs)-1] - NNv[i+(k-1)* (nnew)-1];
        if( (absol <= (thres)) && (absol>= -(thres))  ){
          counti[j-1] = 1;
          meancount = meancount + 1;
        } else {
          counti[j-1] = 0;
        }
      }
      if (meancount >= 1){
        if (normalise >= 1){
          for (j=1; j<=(nobs); ++j){
            WEv[j+(i-1)* (nobs)-1] = WEv[j+(i-1)* (nobs)-1] + counti[j-1]/ (double) (meancount)  ;
          }
        } else {
          for (j=1; j<=(nobs); ++j){
            WEv[j+(i-1)* (nobs)-1] = WEv[j+(i-1)* (nobs)-1] + counti[j-1] ;
          }
        }
      }
    }
  }    
  return WEv ;
} 
