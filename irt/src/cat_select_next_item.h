#ifndef CATSELECTNEXTITEM_H
#define CATSELECTNEXTITEM_H

#include <Rcpp.h>
Rcpp::S4 get_remaining_items(Rcpp::S4 ip, Rcpp::List eh);
double loglik_est_history(Rcpp::List est_history, double theta,
                          bool calculate_loglik = true);
Rcpp::S4 select_next_item_fisher_max_info_cpp(double theta, Rcpp::S4 ip,
                                              int randomesqueN);
Rcpp::List select_next_item_cpp(Rcpp::List cd, Rcpp::List est_history,
                                Rcpp::List additional_args);

#endif
