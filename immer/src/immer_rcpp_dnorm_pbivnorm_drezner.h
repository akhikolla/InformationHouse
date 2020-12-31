//// File Name: immer_rcpp_dnorm_pbivnorm_drezner.h
//// File Version: 1.07

#ifndef _IMMER_IMMER_RCPP_DNORM_PBIVNORM_DREZNER_H
#define _IMMER_IMMER_RCPP_DNORM_PBIVNORM_DREZNER_H
 
#include <Rcpp.h>
using namespace Rcpp;

double const pi0 = 3.14159265359;


double immer_pnorm(double x);

double immer_pnorm(double x);

double immer_signum(double x);

double immer_dnorm2(double a, double b, double rho);

double immer_dnorm(double a);

double pbivnorm_drezner_fct_f_cpp(double x, double y, double ap, double bp, double rho);

double pbivnorm_drezner_all_negative_cpp(double a, double b, double rho);

double pbivnorm_drezner_product_negative_cpp(double a, double b, double rho);

double pbivnorm_drezner_product_positive_cpp(double a, double b, double rho);

double pbivnorm_drezner_numeric_arguments(double a, double b, double rho);

Rcpp::NumericVector pbivnorm_drezner(Rcpp::NumericVector a, Rcpp::NumericVector b,
        Rcpp::NumericVector rho);

double pbivnorm_drezner_derivative_rho_numeric(double a, double b, double rho);

double pbivnorm_drezner_derivative_a_numeric(double a, double b, double rho);

double pbivnorm_drezner_derivative_b_numeric(double a, double b, double rho);

Rcpp::List pbivnorm_drezner_derivative(Rcpp::NumericVector a,
    Rcpp::NumericVector b, Rcpp::NumericVector rho);

double immer_logdnorm2( double x, double y, double mu1, double mu2, double var1,
        double var2, double cov12 );

Rcpp::List immer_logdnorm2_derivative( double x, double y, double mu1, double mu2, double var1,
        double var2, double cov12 );

#endif // _IMMER_IMMER_RCPP_DNORM_PBIVNORM_DREZNER_H
