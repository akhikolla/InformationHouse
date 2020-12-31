#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pairup(NumericVector vec, String type = "less") {
  int n = vec.size();
    
  int index = 0;
  
  if(type == "leq") {
    NumericMatrix m((pow((double)n, 2.0) + n) / 2.0, 2);
    
    for(int i = 0; i < n; i++) {
      for(int j = i; j < n; j++) {
        m(index, 0) = vec[i];
        m(index, 1) = vec[j];
          
        index++;
      }
    }
    
    return(m);
  }
  else if(type == "neq") {
    NumericMatrix m(pow((double)n, 2.0), 2);
    
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        m(index, 0) = vec[i];
        m(index, 1) = vec[j];
          
        index++;
      }
    }
    
    return(m);
  }
  else {
    // Default is 'less'
    NumericMatrix m(n * (n-1) / 2.0, 2);
      
    for(int i = 0; i < n - 1; i++) {
      for(int j = i + 1; j < n; j++) {
        m(index, 0) = vec[i];
        m(index, 1) = vec[j];
            
        index++;
      }
    }
      
    return(m);
  }
}
