#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix bmat_schC(double v1, double v2, double v3, NumericVector section_counts) {
  int students = sum(section_counts);
  
  NumericMatrix bmat(students, students);
  
  // Initialize the matrix with v3
  for(int i = 0; i < bmat.nrow(); i++) {
    for(int j = 0; j < bmat.ncol(); j++) {
      bmat(i, j) = v3;
    }
  }
  
  int index = 0;
  for(int i = 0; i < section_counts.size(); i++) {
    for(int j = 0; j < section_counts[i]; j++) {
      for(int k = 0; k < section_counts[i]; k++) {
        if(j == k) {
          bmat(index + j, index + k) = v1;
        }
        else {
          bmat(index + j, index + k) = v2;
        }
      }
    }
    
    index += section_counts[i];
  }
  
  return bmat;
}

// [[Rcpp::export]]
NumericVector bmat2C(double v1, double v2, NumericVector sec) {
  int students = sum(sec);
  
  NumericMatrix bmat(students, students);
  
  // Initialize the matrix
  for(int i = 0; i < bmat.nrow(); i++) {
    for(int j = 0; j < bmat.ncol(); j++) {
      bmat(i, j) = 0;
    }
  }
  
  int index = 0;
  for(int i = 0; i < sec.size(); i++) {
    
    for(int j = 0; j < sec[i]; j++) {
      for(int k = 0; k < sec[i]; k++) {
        if(j == k) {
          bmat(index + j, index + k) = v1;
        }
        else {
          bmat(index + j, index + k) = v2;
        }
      }
    }
    
    index += sec[i];
  }
  
  return bmat;
}
