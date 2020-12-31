#ifndef RESP_LIK_H
#define RESP_LIK_H

#include <Rcpp.h>
// Respond Likelihood
double resp_lik_bare_item_cpp(double resp, double theta, Rcpp::S4 item);
double resp_lik_bare_testlet_cpp(Rcpp::NumericVector resp, double theta,
                                 Rcpp::S4 testlet);
double resp_lik_bare_itempool_cpp(Rcpp::NumericVector resp, double theta,
                                   Rcpp::S4 ip);
#endif
