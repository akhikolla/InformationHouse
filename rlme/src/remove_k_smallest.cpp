#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector remove_k_smallest(NumericVector v, int k) {
  int *indices =(int *) malloc(sizeof(int) * k);
  
  for(int i = 0; i < k; i++) {
    double smallest = 0;
    int smallest_index = 0;
    
    bool first = true;
    
    for(int j = 0; j < v.size(); j++) {
      // Skip indices that have already been removed
      bool removed = false;
      for(int l = 0; l < i; l++) {
        if(indices[l] == j) {
          removed = true;
          break;
        }
      }
      
      if(removed) continue;
      
      if(first || v[j] < smallest) {
        smallest = v[j];
        smallest_index = j;
      }
      
      first = false;
    }
    
    indices[i] = smallest_index;
  }
  
  NumericVector newv(v.size() - k);
  int newk_index = 0;
  for(int i = 0; i < v.size(); i++) {
    
    bool removed = false;
    for(int l = 0; l < k; l++) {
      if(i == indices[l]) {
        removed = true;
        break;
      }
    }
    
    if(removed) continue;
    
    newv[newk_index] = v[i];
    
    newk_index++;
  }
  
  free(indices);
  
  return(newv);
}
