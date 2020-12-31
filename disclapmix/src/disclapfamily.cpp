#include <Rcpp.h>

using namespace Rcpp;
    
// [[Rcpp::export]]
NumericVector disclapglm_linkfun(NumericVector mu) {
  // returns eta = theta = g(mu) = log(p), where this is g(mu)
  return log( (sqrt(1.0 + mu*mu) - 1.0) / mu);
}

// [[Rcpp::export]]
NumericVector disclapglm_linkinv(NumericVector eta) {
  // returns mu = g^-1(eta), where this is g^-1(eta) = 2*e^eta / (1 - e^(2*eta))
  double eps = std::numeric_limits<double>::epsilon();
  NumericVector expeta = exp(eta);
  return pmax(2*expeta / (1.0 - expeta*expeta), eps);
}

// [[Rcpp::export]]
NumericVector disclapglm_mu_eta(NumericVector eta) {
  double eps = std::numeric_limits<double>::epsilon();
  
  NumericVector exp_eta(exp(eta));
  NumericVector exp_eta_sq(exp_eta*exp_eta);
  NumericVector exp_eta_sq_one(exp_eta_sq - 1.0);
    
  return pmax( (2*exp_eta*(1.0 + exp_eta_sq)) / (exp_eta_sq_one*exp_eta_sq_one), eps);    
}

// [[Rcpp::export]]
NumericVector disclapglm_varfunc(NumericVector mu) {
  return mu * sqrt(1.0 + mu*mu);
}

// [[Rcpp::export]]
double disclapglm_loglikeh(double mu, double y) {
  double p = (mu < 1e-4) ? 0.5*mu : (sqrt(1.0 + mu*mu) - 1.0) / mu;
  double logl = log(1.0-p) - log(1.0+p) + y*log(p);  
  return logl;
}

// [[Rcpp::export]]
double disclapglm_deviance(NumericVector y, NumericVector mu, NumericVector wt) {
  int n = y.size();
  NumericVector dev(n);

  for (int i = 0; i < n; i++) {
    double y_ele = y[i];
    double mu_ele = mu[i];
    
    if (!R_FINITE(mu_ele) || mu_ele < 1e-14) { // Includes NA/NaN
      // Happens if mu is very small
      mu_ele = 1e-14;
    }
    
    if ((int)y_ele == 0) {  
      // y == 0: dev = 2*log((1+p)/(1-p))
      double p = (mu_ele < 1e-4) ? 0.5*mu_ele : (sqrt(1.0 + mu_ele*mu_ele) - 1.0) / mu_ele;
      dev[i] = 2 * log((1.0+p) / (1.0-p));
    } else {
      // y != 0
      dev[i] = 2 * (disclapglm_loglikeh(y_ele, y_ele) - disclapglm_loglikeh(mu_ele, y_ele));
    }
  }  
  
  dev = dev * wt;
  double deviance = std::accumulate(dev.begin(), dev.end(), 0.0);

  return deviance;
}

