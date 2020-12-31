#ifndef PROB_H
#define PROB_H

#include <Rcpp.h>
// prob
double prob_4pm_bare_cpp(double theta, Rcpp::S4 item, int derivative = 0);
Rcpp::NumericVector prob_grm_bare_cpp(double theta, Rcpp::S4 item,
                                      int derivative = 0);
Rcpp::NumericVector prob_gpcm_bare_cpp(double theta, Rcpp::S4 ip,
                                       int derivative = 0);
Rcpp::NumericVector prob_poly_bare_cpp(double theta, Rcpp::S4 item,
                                       int derivative = 0, bool
                                       expected_value = false);
#endif
