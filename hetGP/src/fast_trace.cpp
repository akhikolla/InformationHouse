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

// Compute trace(A * B)
// [[Rcpp::export]]
double fast_trace(NumericMatrix A, NumericMatrix B) {
  double tmp = 0;
  int nrA = A.nrow();
  int ncA = A.ncol();
  const double* ptrA = (const double*) &A(0,0);
  const double* ptrB = (const double*) &B(0,0);
  for(int i = 0; i < nrA; i++){
    ptrA = &A(i,0);
    for(int j = 0; j < ncA; j++, ptrA+=nrA, ptrB++){
      // tmp += A(i,j)*B(j,i);
      tmp += *ptrA * *ptrB;
    }
  }
  return(tmp);
}

// Compute trace(A * B) when A is symmetric
// [[Rcpp::export]]
double trace_sym(NumericMatrix A, NumericMatrix B) {
  double tmp = 0;
  int nrA = A.nrow();
  const double* ptrA = (const double*) &A(0,0);
  const double* ptrB = (const double*) &B(0,0);
  for(int i = 0; i < nrA; i++){
    for(int j = 0; j < nrA; j++, ptrA++, ptrB++){
      tmp += *ptrA * *ptrB;
    }
  }
  return(tmp);
}

// // Compute trace(A * (g t(g))
// // [[Rcpp::export]]
// double fast_trace_g(NumericMatrix A, NumericVector g){
//   double tmp = 0;
//   for(int i = 0; i < A.nrow(); i++){
//     for(int j = 0; j < A.ncol(); j++){
//       tmp += A(i,j) * g(i) * g(j);
//     }
//   }
//   return(tmp);
// }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// /*** R
// A <- matrix(rnorm(10000),100,100)
// B <- matrix(rnorm(10000),100,100)
// print(c(fast_trace(A, B), sum(diag(A %*% B))))
//
//
// library(microbenchmark)
// microbenchmark(fast_trace(A, B), sum(diag(A %*% B)))
// */
