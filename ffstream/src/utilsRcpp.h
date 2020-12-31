#ifndef GUARD_utilsRcpp_h
#define GUARD_utilsRcpp_h

#include "utils.h"
#include<Rcpp.h>

Rcpp::NumericVector computeOneSidedPvalueRcpp(Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector);

Rcpp::NumericVector combineTwoOneSidedPvaluesRcpp(Rcpp::NumericVector, Rcpp::NumericVector);

Rcpp::NumericVector convertPvaluetoCorrectSide(Rcpp::NumericVector);

Rcpp::NumericVector computeTwoSidedPvalueRcpp(Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector);


Rcpp::NumericVector computeStdNormCdf(Rcpp::NumericVector);

Rcpp::NumericVector makeTwoSidedPvalueOneSidedR(Rcpp::NumericVector);

#endif
