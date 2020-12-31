#ifndef CATTERMINATECAT_H
#define CATTERMINATECAT_H

#include <Rcpp.h>
Rcpp::List next_step_cat_cpp(Rcpp::List true_ability,
                             Rcpp::List cd, 
                             Rcpp::Nullable<Rcpp::List> est_history,
                             Rcpp::Nullable<Rcpp::List> additional_args);
Rcpp::List generate_cat_resp_cpp(Rcpp::List true_ability, 
                                 Rcpp::List cd, 
                                 Rcpp::List est_history,
                                 Rcpp::List additional_args);
Rcpp::List est_ability_cat_cpp(Rcpp::List true_ability, 
                               Rcpp::List cd, 
                               Rcpp::List est_history,
                               Rcpp::List additional_args,
                               bool last_estimate = false);
bool terminate_cat_cpp(Rcpp::List true_ability, 
                       Rcpp::List cd, 
                       Rcpp::List est_history,
                       Rcpp::List additional_args);
#endif
