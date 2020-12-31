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
double moy(NumericVector x) {
  double total = 0;
  int n = x.size();
  for (int i=0; i<n; i++) {
    total += x[i];
  }
  return total/n;
}


// [[Rcpp::export]]
NumericVector vestimates(NumericVector x, int binf, int bsup) {
  int n = x.size();
  NumericVector V(bsup-binf+1);
  for (int i=binf; i<bsup+1; i++) {
    NumericVector A(n);
    for (int j=1; j<i+1; j++){
      A[j-1] = log(n*(x[j+i-1]-x[0])/(2*i));
      A[n-i+j-1] =log(n*(x[n-1]-x[n-i+j-i-1])/(2*i));
    }
    for (int k=i+1; k<n-i+1; k++){
      A[k-1] = log(n*(x[k+i-1]-x[k-i-1])/(2*i));
    }
    V[i-binf]=moy(A);
  }
  return V;
}


