#ifndef INFO_H
#define INFO_H

#include <Rcpp.h>
// info
double info_item_bare_cpp(double theta, Rcpp::S4 item, bool observed, double resp);
Rcpp::NumericVector info_itempool_bare_cpp(double theta, Rcpp::S4 ip, bool tif,
  bool observed, Rcpp::Nullable<Rcpp::NumericVector> resp = R_NilValue);
Rcpp::NumericMatrix info_itempool_cpp(Rcpp::NumericVector theta, Rcpp::S4 ip, bool tif,
  bool observed, Rcpp::Nullable<Rcpp::NumericMatrix> resp);

#endif
