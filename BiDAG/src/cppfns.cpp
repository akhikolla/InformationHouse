#include <Rcpp.h>
using namespace Rcpp;
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
NumericVector extractT(IntegerVector xs, IntegerVector ys, NumericMatrix ts) {
  int ln = ys.size();
  NumericVector out(ln);
  
  for(int i = 0; i < ln; ++i) {
    out[i] = ts(xs[i],ys[i]);
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix collectCcatwt(IntegerVector xs, IntegerVector ys, NumericVector ws, int n, int m) {
  int ln = ys.size();
  NumericMatrix out(n, m);
  
  for(int i = 0; i < ln; ++i) {
    out(xs[i], ys[i]) += ws[i];
  }
  return out;
}

// [[Rcpp::export]]
IntegerMatrix collectCcat(IntegerVector xs, IntegerVector ys, int n, int m) {
  int ln = ys.size();
  IntegerMatrix out(n, m);
  
  for(int i = 0; i < ln; ++i) {
    out(xs[i], ys[i]) ++;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector collectC(IntegerVector xs, NumericVector ys, int n) {
  int ln = ys.size();
  NumericVector out(n);
  
  for(int i = 0; i < ln; ++i) {
    out[xs[i]] += ys[i];
  }
  return out;
}



// [[Rcpp::export]]
NumericVector takefirst(NumericVector xs, int pos) {
  NumericVector out(pos);
  
  for(int i = 0; i < pos; ++i) {
    out[i] = xs[i];
  }
  return out;
}

// [[Rcpp::export]]
NumericVector takelast(NumericVector xs, int pos, int n) {
  
  NumericVector out(n-pos);
  
  for(int i = pos; i < n; ++i) {
    out[i-pos] = xs[i];
  }
  return out;
}



