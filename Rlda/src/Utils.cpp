#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "progress.hpp"
#include <iostream>
#include <ctime>
#include <fstream>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace arma;
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/


// [[Rcpp::export]]
NumericMatrix convertSBtoNormal(NumericMatrix vmat,
                                int ncol, int nrow,
                                NumericVector prod) {
  NumericMatrix res(nrow,ncol);

  for(int j=0; j<ncol;j++){
    res(_,j)=vmat(_,j)*prod;
    prod=prod*(1-vmat(_,j));
  }

  return (res);
}

// [[Rcpp::export]]
NumericVector aggregatesum(NumericVector Tobesum,
                           int nind, int nobs,
                           IntegerVector ind) {
  NumericVector res(nind);

  for(int i=0; i<nobs; i++){
    for (int j=0; j<nind; j++){
      if (ind[i]==j+1) res[j]=res[j]+Tobesum[i];
    }
  }

  return (res);
}


