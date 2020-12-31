#include <RcppArmadillo.h>

//This code is modified from the repository: https://github.com/coatless/r-to-armadillo which is a collection of R functions written in C++

//[[Rcpp::export]]
arma::vec cfilter(arma::vec x, arma::vec filter)
{
  
  int nx = x.n_elem;
  int nf = filter.n_elem;
  int nshift = nf/2;
  
  double z, tmp;
  
  
  arma::vec out = arma::zeros<arma::vec>(nx);
  
  for(int i = 0; i < nx; i++) {
    z = 0;
    if(i + nshift - (nf - 1) < 0 || i + nshift >= nx) {
      out(i) = NA_REAL;
      continue;
    }
    for(int j = std::max(0, nshift + i - nx); j < std::min(nf, i + nshift + 1) ; j++) {
      tmp = x(i + nshift - j);
      z += filter(j) * tmp;
    }
    out(i) = z;
  }
  
  return out;
}


//[[Rcpp::export]]
arma::vec mldivide(arma::mat A, arma::vec B){
  return(arma::solve(A,B));
}
