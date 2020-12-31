// functions that remember against which special external libraries
// the package transport has been linked

#include <Rcpp.h>

using namespace Rcpp;


//[[Rcpp::export]]
SEXP openmp_present(){
#ifdef _OPENMP
  return Rcpp::wrap(1);
#endif
#ifndef _OPENMP
  return Rcpp::wrap(0);
#endif
}


// [[Rcpp::export]]
SEXP cplex_present() {
#ifdef WITH_CPLEX
  return Rcpp::wrap(1);
#endif
#ifndef WITH_CPLEX
  return Rcpp::wrap(0);
#endif
}


// [[Rcpp::export]]
SEXP cgal_present() {
#ifdef WITH_CGAL
  return wrap(1);
#endif
#ifndef WITH_CGAL
  return wrap(0);
#endif
}
