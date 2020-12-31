#include <Rcpp.h>
using namespace Rcpp;

//' C++ Helper for \code{episode_check()}
//'
//' Meant to be a C++ helper within episode_check() to improve speed
//
//' @param x Our gap vector
//' @param perm_gap An integer value representing the permissible gap in therapy, if the gap is greater than this number, a new episode will be created
//' @param init_rank An integer representing the initial ranking (What should we label episode 1?). The default is 1. 
//' @return A new vector labeling each rows corresponding episode of care to be appended to a users dataframe


// [[Rcpp::export(name = ".episode_checkCpp")]]
NumericVector episode_checkCpp(NumericVector x, int perm_gap, int init_rank){
  int data_length = x.size();
  NumericVector out(data_length);
  
  for(int i = 0; i < data_length; ++i){
    
    // if the permissible gap is less than gap... 
    if(perm_gap < x[i]){
      // update the current ranking
      init_rank = init_rank + 1;
    }
    // append current ranking to out
    out[i] = init_rank;
  }
  return out;
}

