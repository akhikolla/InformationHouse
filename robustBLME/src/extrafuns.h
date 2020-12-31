#ifndef EXTRAFUNS_HPP_
#define EXTRAFUNS_HPP_

#include <RcppArmadillo.h>
using namespace Rcpp;


double psi_huber(double x, const double c);

arma::vec vpsi_huber(arma::vec x, const double c, int xLen);

double psip_huber(double x, const double c);

double dhalfCauchy(double x, double scale, bool lg = false);

double
  dinvgamma (double x, double shape, double scale, bool lg);

arma::mat thinMat(arma::mat X, arma::vec index);

arma::vec rmvnorm(arma::vec mu, arma::mat S, int p);

arma::vec rmvnorm2(arma::vec mu, arma::mat S, int p);

arma::vec rmvt(arma::vec mu, arma::mat S, int p, double df);

double
  dmvt(arma::vec x, arma::vec mu, arma::mat S, int p, double df, bool lg);

#endif //EXTRAFUNS_HPP_
