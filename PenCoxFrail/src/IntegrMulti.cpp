// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List IntegrMultiCpp(NumericMatrix alpha2,
                    NumericVector linvec2,
              NumericMatrix Phibig2,
              NumericMatrix tUeach2,
              NumericMatrix Ueach2,
              NumericMatrix Tmat2,
              NumericMatrix Phidoublebig2,
              NumericMatrix Udoubleeach2,
              NumericMatrix XW2) { 
  
  mat alpha = as<arma::mat>(alpha2);
  mat Phibig = as<arma::mat>(Phibig2);
  mat tUeach = as<arma::mat>(tUeach2);
  mat Ueach = as<arma::mat>(Ueach2);
  mat Tmat = as<arma::mat>(Tmat2);
  mat Phidoublebig = as<arma::mat>(Phidoublebig2);
  mat Udoubleeach = as<arma::mat>(Udoubleeach2);
  mat XW = as<arma::mat>(XW2);
  colvec linvec = as<arma::colvec>(linvec2);
  
  
  colvec linpred = XW * linvec;

  colvec eta = exp(linpred); 

  rowvec teta = trans(eta);

  mat helpcalc = trans(exp(Phibig * ( tUeach % alpha)) % Tmat);

  mat intma = (helpcalc * Phibig) % Ueach;

  rowvec intarray = teta * ((helpcalc * Phidoublebig) % Udoubleeach);

  return List::create(Named("int.array") = intarray,
                      Named("int.ma") = intma,
                      Named("eta") = eta,
                      Named("lin.pred") = linpred,
                      Named("t.eta") = teta,
                      Named("help.calc") = helpcalc);
}

