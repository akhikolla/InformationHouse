
#include <Rcpp.h>
#include "typedefs.h"
#include "headers.h"

using namespace Rcpp;


// allTRUE and allFALSE
// --------------------

//// [[Rcpp::export]]
bool C_allTRUE(const LogicalVector x){
  int len = x.size();
  bool out = true;
  for (int i=0; i<len; i++){
    if (!x[i]){
      out = false;
      break;
    }
  }
  return out;
}
//// [[Rcpp::export]]
bool C_allFALSE(const LogicalVector x){
  int len = x.size();
  bool out = true;
  for (int i=0; i<len; i++){
    if (x[i]){
      out = false;
      break;
    }
  }
  return out;
}


// Functions to manipulate intList's
// ---------------------------------

// Auxiliary fun used in isSubsetOf 
// IntegerVector rmFirst(const IntegerVector x) {
//   int n = x.size();
//   IntegerVector out(n-1);
//   if (n>1) out = x[Range(1, n-1)];
//   return out;
// }

// [[Rcpp::export]]
bool C_isSubsetOf(const IntegerVector x, const IntegerVector y) {
  int nx = x.size();
  int ny = y.size();
  int i=0, j=0;
  bool out;
  while (true) {
    if (i == nx){  // reached end of in x? => TRUE
      out = true;
      break;
    }
    if (j == ny){  // reached end of in y (but not yet end of x)? => FALSE
      out = false;
      break;
    }
    if (nx-i > ny-j){ // more elements left in x than in y? => FALSE
      out = false;
      break;
    }
    int x1=x[i], y1=y[j];
    if (x1 < y1){     // x_i < y_j => FALSE
      out = false;
      break;
    }
    if (x1 == y1){
      ++i;
    }
    ++j;
  }
  return(out);
}


// Check if two intList's are equivalent
// [[Rcpp::export]]
bool intList_equal(const intList x, const intList y){
  int lenx = x.size(), leny = y.size();
  if (lenx != leny) return(false);
  bool found=false;
  for (int i=0; i<lenx; ++i){
    found=false;
    for (int j=0; j<leny; ++j){
      if (setequal(as<IntegerVector>(x[i]), as<IntegerVector>(y[j]))) found=true;
      if (found) continue;
    }
    if (!found) break;
  }
  return(found);
}

// [[Rcpp::export]]
LogicalVector C_hasSubsetIn(const intList y, const intList x){
  int ylen = y.size(), xlen = x.size();
  LogicalVector out(ylen);
  for (int iy=0; iy < ylen; iy++) {
    bool foundSubset=false;
    for (int ix=0; ix < xlen; ix++) {
      if (C_isSubsetOf(x[ix], y[iy])){
        foundSubset=true;
        break;
      }
    }
    out[iy]=foundSubset;
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector C_hasSupersetIn(const intList x, const intList y, const bool ignore_equals){
  int ylen = y.size(), xlen = x.size();
  LogicalVector out(xlen);
  for (int ix=0; ix<xlen; ix++) {
    bool foundSuperset=false;
    for (int iy=0; iy<ylen; iy++) {
      if (C_isSubsetOf(x[ix], y[iy])){
        if (ignore_equals && setequal(x[ix], y[iy])){
          continue;
        } else {
          foundSuperset=true;
          break;
        }
      }
    }    
    out[ix]=foundSuperset;
  }  
  return out;
}


// Version of C_hasSubsetIn that appies to IntegerMatrix'es
// [[Rcpp::export]]
LogicalVector C_hasSubsetInM(const IntegerMatrix y, const IntegerMatrix x){ 
  // "M" indicates application to matrices
  int nry=y.nrow(), nrx=x.nrow();
  LogicalVector out(nry);
  for (int iy=0; iy < nry; iy++) {
    IntegerVector y1 = y(iy, _);
    bool foundSubset=false;
    for (int ix=0; ix < nrx; ix++) {
      IntegerVector x1 = x(ix, _);
      if (C_isSubsetOf(x1, y1)){
        foundSubset=true;
        break;
      }
    }    
    out[iy]=foundSubset;
  }  
  return(out);
}


/*/ Print an intlist to the console
// [[Rcpp::export]]
void C_show_intList(intList x){
  int lenx = x.size();
  for (int i=0; i<lenx; i++){
    IntegerVector xi = x[i];
    Rcpp::Rcout << xi << " ";
    if (i<lenx-1) Rcpp::Rcout << "/ ";
  }
  Rcpp::Rcout << std::endl;
  return;
}*/


// Append two intList's
// [[Rcpp::export]]
intList C_append_intList(const intList x, const intList y){
  int lenx = x.size(), leny = y.size();
  List out(lenx + leny);
  for (int i=0; i<lenx; i++){
    out[i] = x[i];
  }
  for (int i=0; i<leny; i++){
    out[lenx + i] = y[i];
  }
  return(as<intList>(out));
}

