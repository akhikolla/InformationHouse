#ifndef ABILITYESTIMATION_H
#define ABILITYESTIMATION_H

#include <Rcpp.h>
double est_ability_4pm_nr_itempool_cpp(
  Rcpp::NumericVector resp, Rcpp::S4 ip, Rcpp::NumericVector theta_range, 
  double criterion = 0.001, 
  Rcpp::Nullable<Rcpp::NumericVector> initial_estimates = R_NilValue);
Rcpp::List est_ability_eap_single_examinee_cpp(
    Rcpp::NumericVector resp, Rcpp::S4 ip, Rcpp::NumericVector theta_range,
    int no_of_quadrature, Rcpp::String prior_dist, Rcpp::NumericVector prior_par);
Rcpp::List est_ability_owen_cpp(Rcpp::S4 ip, Rcpp::NumericVector resp,
                                double m0, double v0);
#endif
