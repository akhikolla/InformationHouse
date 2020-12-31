#ifndef RESPLOGLIK_H
#define RESPLOGLIK_H

#include <Rcpp.h>
// Respond Log-Likelihood
double resp_loglik_bare_item_cpp(double resp, double theta, Rcpp::S4 item, 
                                 int derivative = 0);
double resp_loglik_bare_testlet_cpp(Rcpp::NumericVector resp, double theta,
                                    Rcpp::S4 testlet, int derivative = 0);
double resp_loglik_bare_itempool_cpp(Rcpp::NumericVector resp, double theta,
                                      Rcpp::S4 ip, int derivative = 0);
#endif
