#ifndef GUARD_utilsRcpp_cpp
#define GUARD_utilsRcpp_cpp

#include "utils.h"
#include<Rcpp.h>

//this class just exposes functions in utils.cpp to Rcpp


// [[Rcpp::export]]
Rcpp::NumericVector computeOneSidedPvalueRcpp(Rcpp::NumericVector x_, Rcpp::NumericVector a_, Rcpp::NumericVector b_){
    double x = Rcpp::as<double> (x_);
    double a = Rcpp::as<double> (a_);
    double b = Rcpp::as<double> (b_);
    double ans = computeOneSidedPvalue(x, a, b);

    return Rcpp::NumericVector::create(ans);
}


// [[Rcpp::export]]
Rcpp::NumericVector combineTwoOneSidedPvaluesRcpp(Rcpp::NumericVector p1_, Rcpp::NumericVector p2_){
    double p1 = Rcpp::as<double> (p1_);
    double p2 = Rcpp::as<double> (p2_);
    double ans = combineTwoOneSidedPvalues(p1, p2);

    return Rcpp::NumericVector::create(ans);
}




// [[Rcpp::export]]
Rcpp::NumericVector convertPvalueToCorrectSideRcpp(Rcpp::NumericVector p_){
    double p = Rcpp::as<double> (p_);
    double ans = convertPvalueToCorrectSide(p);

    return Rcpp::NumericVector::create(ans);
}




// [[Rcpp::export]]
Rcpp::NumericVector computeTwoSidedPvalueRcpp(Rcpp::NumericVector x_, Rcpp::NumericVector a_, Rcpp::NumericVector b_){
    double x = Rcpp::as<double> (x_);
    double a = Rcpp::as<double> (a_);
    double b = Rcpp::as<double> (b_);
    double ans = computeTwoSidedPvalue(x, a, b);

    return Rcpp::NumericVector::create(ans);
}


// [[Rcpp::export]]
Rcpp::NumericVector computeStdNormCdf(Rcpp::NumericVector x_){
    double x = Rcpp::as<double> (x_);
    double ans = stdnormcdf(x);

    return Rcpp::NumericVector::create(ans);
}


// [[Rcpp::export]]
Rcpp::NumericVector makeTwoSidedPvalueOneSidedR(Rcpp::NumericVector p2_){
    double p2 = Rcpp::as<double> (p2_);
    return Rcpp::NumericVector::create(makeTwoSidedPvalueOneSided(p2));
}

#endif
