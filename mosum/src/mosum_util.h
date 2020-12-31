#ifndef __MOSUM_UTIL_H
#define __MOSUM_UTIL_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector rolling_sum(const NumericVector &x, unsigned G);

IntegerVector eta_criterion_help(const IntegerVector &candidates, 
                                 const NumericVector &m_values,
                                 double eta, double G_left, double G_right);

#endif
