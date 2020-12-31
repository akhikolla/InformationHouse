#include <Rcpp.h>
#include <string>
#include "typedefs.h"

using namespace Rcpp;

// Some functions for manipulating character vectors
// -------------------------------------------------

std::string concat2(const std::string x, const std::string y, const std::string sep){
  return x + sep + y;
}
// paste(x, collapse = sep)
// [[Rcpp::export]]
std::string C_concat(const CharacterVector x, const std::string sep){
  std::string out = ""; 
  std::string sep0 = "";
  for (int i=0; i<x.length(); i++){
    if (i==1){
      sep0 = sep;
    }
    out = concat2(out, Rcpp::as<std::string>(x[i]), sep0);
  }
  return out;
}

// applies C_concat to each element of x separately 
// optionally sorts the elements (disjuncts)
// [[Rcpp::export]]
CharacterVector C_mconcat(const charList x, const std::string sep, const bool sorted = false){
  int n = x.size();
  CharacterVector out(n);
  for (int i=0; i<n; i++){
    CharacterVector xi = x[i];
    CharacterVector xis = Rcpp::clone(xi);
    if (sorted){
      xi = xis.sort();
    }
    out(i) = C_concat(xis, sep);
  }
  return out;
}

// convert a charlist to a character string
// [[Rcpp::export]]
std::string C_charList2string(const charList x, const std::string disj = "+", 
                              const std::string conj = "*", const bool sorted = false){
  CharacterVector conjs = C_mconcat(x, conj, sorted);
  if (sorted){
    conjs = conjs.sort();
  }
  return C_concat(conjs, disj);
}

// convert a recursive charlist to a character vector
// [[Rcpp::export]]
CharacterVector C_recCharList2char(const recCharList x, const std::string disj = "+", 
                                   const std::string conj = "*", const bool sorted = false){
  int n = x.size();
  CharacterVector out(n);
  for (int i=0; i<n; i++){
    charList xi = as<charList>(x[i]);
    out[i] = C_charList2string(xi, disj, conj, sorted);
  }
  return out;
}

