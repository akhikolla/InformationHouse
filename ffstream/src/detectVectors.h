#ifndef GUARD_detectvectors_h
#define GUARD_detectvectors_h

#include "detector.h"
#include "fff.h"
#include "aff.h"
#include "fffcd.h"
#include "affcd.h"
#include "cusumcd.h"
#include "ewmacd.h"
//#include<Rcpp.h>
// #include<vector>


//FFF methods
Rcpp::List cpp_detectFFFMeanMultiple(Rcpp::NumericVector, double, double, int);
Rcpp::List cpp_detectFFFMeanSingle(Rcpp::NumericVector, double, double, int);
Rcpp::List cpp_detectFFFMeanSinglePrechange(Rcpp::NumericVector, double, double, double, double);

//AFF methods
Rcpp::List cpp_detectAFFMeanMultiple(Rcpp::NumericVector, double, double, int);
Rcpp::List cpp_detectAFFMeanSingle(Rcpp::NumericVector, double, double, int);
Rcpp::List cpp_detectAFFMeanSinglePrechange(Rcpp::NumericVector, double, double, double, double);

//CUSUM methods
Rcpp::List cpp_detectCUSUMMeanMultiple(Rcpp::NumericVector, double, double, int);
Rcpp::List cpp_detectCUSUMMeanSingle(Rcpp::NumericVector, double, double, int);
Rcpp::List cpp_detectCUSUMMeanSinglePrechange(Rcpp::NumericVector, double, double, double, double);

//EWMA methods
Rcpp::List cpp_detectEWMAMeanMultiple(Rcpp::NumericVector, double, double, int);
Rcpp::List cpp_detectEWMAMeanSingle(Rcpp::NumericVector, double, double, int);
Rcpp::List cpp_detectEWMAMeanSinglePrechange(Rcpp::NumericVector, double, double, double, double);


// double cpp_detectAFFMean(Rcpp::NumericVector, double);


#endif
