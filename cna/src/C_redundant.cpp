
#include <Rcpp.h>
#include "headers.h"

using namespace Rcpp;

// min(which(!x))-1
// returns 0 if all(x) is TRUE
// [[Rcpp::export]]
int C_find_first_false(const LogicalVector x){
  int i = 0;
  while(x(i)){
    i++;
    if (i == x.size()){
      i = -1;
      break;
    }
  }
  return i;
}

// selects the rows with only one FALSE and determines the columns with only TRUE value in these rows
// [[Rcpp::export]]
LogicalVector C_redund(const LogicalMatrix x){
  int k = x.ncol();
  LogicalVector possibly_redundant = rep(true, k);
  for (int i=0; i < x.nrow(); i++){
    LogicalVector xi = x(i, _);
    if (sum(as<IntegerVector>(xi)) != k-1) continue;
    int pos = C_find_first_false(xi);
    possibly_redundant(pos) = false;
    if (C_allFALSE(possibly_redundant)) break;
  }
  return possibly_redundant;
}


// Applies C_redund repeatedly
// l specifies the grouping of the columns in x - thus sum(l) must equal ncol(x)
// [[Rcpp::export]]
List C_mredund(const LogicalMatrix x, const IntegerVector l){
  if (sum(l) != x.ncol()) stop("sum(l) must equal ncol(x)");
  int start, end=-1;
  List out(l.size());
  for (int i=0; i<l.size(); i++){
    start = end+1;
    end += l[i]; 
    LogicalMatrix subx = clone(x)(_, Range(start, end));
    out[i] = C_redund(subx);
  }
  return out;
}
