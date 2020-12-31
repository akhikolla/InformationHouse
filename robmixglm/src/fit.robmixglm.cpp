#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector llrandpoiscpp(NumericVector y, NumericVector lp, NumericVector tau2, NumericMatrix gh) {

int nquad=gh.nrow();
int nobs=y.size();

NumericVector ll(nobs);
NumericVector thevals(nquad);
double thel, maxval;

for (int irow=0;irow<nobs;irow++) {
  if (tau2[0]==0.0) ll(irow) = R::dpois(y(irow), std::exp(lp(irow)), true);
  else {
    for (int j=0;j<nquad;j++) {
       thevals(j)=R::dpois(y(irow),std::exp(lp(irow)+gh(j,0)*sqrt(tau2[0])),true)+std::log(gh(j,1));
    }
    maxval=max(thevals);
    thel=0.0;
    for (int j=0;j<nquad;j++) {
      thel=thel+std::exp(thevals(j)-maxval);
    }
    ll(irow) = maxval+std::log(thel);
  }
}
return ll;
}


// [[Rcpp::export]]
NumericVector llrandtruncpoiscpp(NumericVector y, NumericVector lp, NumericVector tau2, NumericMatrix gh) {
    
    int nquad=gh.nrow();
int nobs=y.size();

NumericVector ll(nobs);
NumericVector thevals(nquad);
double thel, maxval;

for (int irow=0;irow<nobs;irow++) {
  /* ???? call actuar::dztpois
  */
  if (tau2[0]==0.0) ll(irow) = R::dpois(y(irow), std::exp(lp(irow)), true)-std::log(1-std::exp(-(std::exp(lp(irow)))));
  else {
    for (int j=0;j<nquad;j++) {
      double lambda = std::exp(lp(irow)+gh(j,0)*sqrt(tau2[0]));
      thevals(j)=R::dpois(y(irow), lambda, true)-std::log(1-std::exp(-lambda))+std::log(gh(j,1));
    }
    maxval=max(thevals);
    thel=0.0;
    for (int j=0;j<nquad;j++) {
      thel=thel+std::exp(thevals(j)-maxval);
    }
    ll(irow) = maxval+std::log(thel);
  }
}
return ll;
}

// [[Rcpp::export]]
NumericVector llrandgammacpp(NumericVector y, NumericVector lp, NumericVector tau2, NumericVector phi, NumericMatrix gh) {

int nquad=gh.nrow();
int nobs=y.size();

NumericVector ll(nobs);
NumericVector thevals(nquad);
double thel, maxval;

for (int irow=0;irow<nobs;irow++) {
  thel=0.0;
  if (tau2[0]==0.0) ll(irow) = R::dgamma(y(irow), 1.0/phi[0], 1.0/(phi[0]*std::exp(lp(irow))), true);
  else {
    for (int j=0;j<nquad;j++) {
      thevals(j)=R::dgamma(y(irow), 1.0/phi[0], 1.0*(phi[0]*std::exp(lp(irow)+gh(j,0)*sqrt(tau2[0]))),true)+std::log(gh(j,1));
    }
    maxval=max(thevals);
    thel=0.0;
    for (int j=0;j<nquad;j++) {
      thel=thel+std::exp(thevals(j)-maxval);
    }
    ll(irow) = maxval+std::log(thel);
  }
}
return ll;
}

// [[Rcpp::export]]
NumericVector llrandbinomcpp(NumericMatrix y, NumericVector lp, NumericVector tau2, NumericMatrix gh) {

int nquad=gh.nrow();
int nobs=y.nrow();

NumericVector ll(nobs);
NumericVector thevals(nquad);
double thel, maxval;

for (int irow=0;irow<nobs;irow++) {
  thel=0.0;
  if (tau2[0]==0.0) ll(irow) = R::dbinom(y(irow,0),y(irow,0)+y(irow,1), 1.0/(1.0+std::exp(-lp(irow))), true);
  else {
    for (int j=0;j<nquad;j++) {
      thevals(j)=R::dbinom(y(irow,0),y(irow,0)+y(irow,1),1.0/(1.0+std::exp(-lp(irow)-gh(j,0)*sqrt(tau2[0]))),true)+std::log(gh(j,1));
    }
    maxval=max(thevals);
    thel=0.0;
    for (int j=0;j<nquad;j++) {
      thel=thel+std::exp(thevals(j)-maxval);
    }
    ll(irow) = maxval+std::log(thel);
  }
}
return ll;
}
