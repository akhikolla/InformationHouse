
#include <Rcpp.h>
#include "typedefs.h"
#include "headers.h"

using namespace Rcpp;

//=============================================================================

// C++-Functions related to is.minimal and is.submodel
// ---------------------------------------------------

// Determine if an intList x is minimal with respect to a reference intList
// (ref contains existing solutions)
// if shortcut is true, (row-wise) execution is stopped as soon as a complete line
// of false values is encountered
// [[Rcpp::export]]
LogicalMatrix C_disj_contained(const intList x, const intList y, const bool shortcut){
  int lenx = x.size(), leny = y.size();
  LogicalMatrix out(lenx, leny);
  for (int i=0; i<lenx; ++i){
    bool anyTrueInRow=false;
    for (int j=0; j<leny; ++j){
      out(i,j) = C_isSubsetOf(as<IntegerVector>(x[i]), as<IntegerVector>(y[j]));
      // Rcout << i << " " << j << ": " << out(i,j) << std::endl;
      if (shortcut & !anyTrueInRow & out(i,j)) anyTrueInRow = true;
    }
    if (shortcut & !anyTrueInRow) break;  // if complete row without TRUE -> break
  }
  return(out);
}
 

// [[Rcpp::export]]
LogicalVector C_is_submodel(const recIntList x, const intList ref, const bool strict){
  int lenx = x.size();
  LogicalVector out(lenx);
  for (int i=0; i<lenx; ++i){
    out[i] = C_checkHallsCondition(C_disj_contained(as<intList>(x[i]), ref, true));
    if (out[i] & strict & intList_equal(as<intList>(x[i]), ref)){
      out[i]=false;
    }
  }
  return(out);
}


// Determine the intLists in x that are minimal with respect to existing solutions ref
// Use separate functions for cases strictly and not strictly minimal
// [[Rcpp::export]]
LogicalVector C_minimal(const recIntList x, const recIntList ref, 
                        bool strict){ 
  int lenx = x.size(), lenr = ref.size();
  LogicalVector out(lenx);
  for (int i=0; i<lenx; ++i){
    bool found=false;
    for (int j=0; j<lenr; ++j){
      found = C_checkHallsCondition(
                C_disj_contained(as<intList>(ref[j]), as<intList>(x[i]), true));
      if (found & strict & 
          intList_equal(as<intList>(x[i]), as<intList>(ref[j]))){
        found=false;
      }
      if (found) break;
    }
    out[i] = !found;
  }
  return(out);
}

// Functions (auxiliary and main) to check Hall's condition for minimality
// -----------------------------------------------------------------------

//// [[Rcpp::export]]
IntegerVector initComb(const int k){
  return seq(0, k-1);
}
//// [[Rcpp::export]]
void nextComb(IntegerVector ii, const int k, const int n){
  if (k == 1){
    ++ii[0];
    return;
  }
  ++ii[k-1];
  if (ii[k-1] == n){
    nextComb(ii, k-1, n-1);
    ii[k-1] = (ii[k-2]+1) % n;
  }
}
//// [[Rcpp::export]]
bool checkLastComb(const IntegerVector ii, const int k, const int n){
  return (ii[0] == n-k);
}

// colAnys(x[rows, ])
// [[Rcpp::export]]
LogicalVector C_rowSubsetColAnys(const LogicalMatrix x, const IntegerVector rows){
  int nc=x.ncol(), k=rows.size(); 
  LogicalVector out(nc);
  for (int i=0; i<k; ++i){
    int r=rows[i];
    for (int j=0; j<nc; ++j){
      out[j] |= static_cast<bool>(x(r,j));
    }
  }
  return out;
}

// Hall's condition for one size of subset
// [[Rcpp::export]]
bool C_checkHall_k(const LogicalMatrix x, const int k) {
  bool out=true;
  int n=x.nrow();
  if (k>n) stop("k too large");
  IntegerVector ii = initComb(k);
  do {
    LogicalVector colAnys=C_rowSubsetColAnys(x, ii);
    // Rcpp::Rcout << "ii: " << ii << std::endl;
    // Rcpp::Rcout << "colAnys: " << colAnys << std::endl;
    if (sum(as<IntegerVector>(colAnys)) < k){
      out=false;
    } else {
      if (checkLastComb(ii, k, n)) break;
      nextComb(ii, k, n);
    }
  } while (out);
  return out;
}

// Hall's condition overall
// [[Rcpp::export]]
bool C_checkHallsCondition(const LogicalMatrix x){
  int n=x.nrow();
  bool out=true;
  for (int i=0; i<n; ++i){
    if (!C_checkHall_k(x, i+1)){
      out=false;
      break;
    }
  }
  return out;
}


    // Old (incorrect) version of C_minimal
    // ------------------------------------

    // Determine if an intList x is minimal with respect to a reference intList
    // (ref contains existing solutions)
    // [[Rcpp::export]]
    bool C_intList_minimal_old(const intList x, const recIntList ref, bool ignore_equals){
      int lenr = ref.size();
      bool is_minimal = true;
      for (int j=0; j<lenr; j++){
        //Rcpp::Rcout << "j: " << j << std::endl;
        LogicalVector lv = C_hasSupersetIn(as<intList>(ref[j]), x, false);
        //Rcpp::Rcout << "lv: " << lv << std::endl;
        if (C_allTRUE(lv)){ 
          is_minimal = false;
          if (ignore_equals){
            bool eq = C_allTRUE(C_hasSupersetIn(x, as<intList>(ref[j]), false));
            if (eq){ 
              is_minimal = true;
            }
          }
          //Rcpp::Rcout << "is_minimal: " << is_minimal << std::endl;
          if (!is_minimal) break;
        }
      }
      return(is_minimal);
    }
    
    // Determine the intLists in x that are minimal with respect to existing solutions ref
    // [[Rcpp::export]]
    LogicalVector C_minimal_old(const recIntList x, const recIntList ref, bool ignore_equals){
      int lenx = x.size();
      LogicalVector out(lenx);
      for (int i=0; i<lenx; i++){
        out[i] = C_intList_minimal_old(as<intList>(x[i]), ref, ignore_equals);  
      }
      return(out);
    }
