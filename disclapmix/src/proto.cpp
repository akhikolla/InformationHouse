#include <Rcpp.h>

#include "eigenmvn.h"
#ifndef M_PI
#define M_PI REAL(3.1415926535897932384626433832795029)
#endif

using namespace Rcpp;

/*
// [[Rcpp::export]]
NumericVector rcpp_calculate_haplotype_probabilities(IntegerMatrix new_data, IntegerMatrix x, IntegerMatrix y, NumericMatrix p, NumericVector tau, int nsim = 1000) {
  if (nsim <= 0) {
    stop("nsim must be at least 1")
  }
  
  Eigen::VectorXd mean;
  Eigen::MatrixXd covar;
  

  Eigen::EigenMultivariateNormal<double> normX_solver(mean, covar);
  
  MatrixXd = normX_solver.samples(5000);}

*/
