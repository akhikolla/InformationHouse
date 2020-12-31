#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector blocksums_i(IntegerVector x, int window) {
  
  // Get length(x) and initialize output vector
  int n = x.size();
  NumericVector out(n / window);
  
  // Loop through and calculate block sums
  double sum = 0;
  int index = 0;
  int resetindex = 1;
  for (int a = 0; a < n; ++a) {
    sum += x[a];
    if (resetindex == window) {
      out[index] = sum;
      index += 1;
      sum = 0;
      resetindex = 0;
    }
    resetindex += 1;
  }              
  return out;
  
}

// [[Rcpp::export]]
NumericVector blocksums_i_max(IntegerVector x, int window) {
  
  // Get length(x)
  int n = x.size();
  
  // First block sum
  double sum = 0;
  for (int a = 0; a < window; ++a) {
    sum += x[a];
  }
  
  // Loop through and find max block sum
  double max = sum;
  sum = 0;
  int resetindex = 1;
  for (int a = window; a < n; ++a) {
    sum += x[a];
    if (resetindex == window) {
      if (sum > max) max = sum;
      sum = 0;
      resetindex = 0;
    }
    resetindex += 1;
  }
  return max;
  
}

// [[Rcpp::export]]
NumericVector blocksums_n(NumericVector x, int window) {
  
  // Get length(x) and initialize output vector
  int n = x.size();
  NumericVector out(n / window);
  
  // Loop through and calculate block sums
  double sum = 0;
  int index = 0;
  int resetindex = 1;
  for (int a = 0; a < n; ++a) {
    sum += x[a];
    if (resetindex == window) {
      out[index] = sum;
      index += 1;
      sum = 0;
      resetindex = 0;
    }
    resetindex += 1;
  }              
  return out;
  
}

// [[Rcpp::export]]
NumericVector blocksums_n_max(NumericVector x, int window) {
  
  // Get length(x)
  int n = x.size();
  
  // First block sum
  double sum = 0;
  for (int a = 0; a < window; ++a) {
    sum += x[a];
  }
  
  // Loop through and find max block sum
  double max = sum;
  sum = 0;
  int resetindex = 1;
  for (int a = window; a < n; ++a) {
    sum += x[a];
    if (resetindex == window) {
      if (sum > max) max = sum;
      sum = 0;
      resetindex = 0;
    }
    resetindex += 1;
  }
  return max;
  
}
