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
NumericVector repeatgenerate(NumericMatrix W, double lambda) {
  int i,j;
  int nrow=W.nrow();
  int ncol=W.ncol();
  NumericVector mi(nrow),ci(nrow),Wi(nrow),mean_eib(nrow),sd_eib(nrow);
  NumericVector Temp(nrow*ncol);
  NumericVector out(nrow);
  NumericMatrix Tbij(nrow,ncol);
  
  double Ws;
  
  for(i = 0; i < nrow; ++i){
    Ws = 0;
    for (j = 0; j < ncol; ++j){
      if (NumericMatrix::is_na(W(i,j))==FALSE) {
        mi(i) += 1;
        Ws += W(i,j);}
    }
    ci(i) = sqrt(lambda/mi(i)/(mi(i)-1));
    Wi(i) = Ws/mi(i);
  }
  Temp = rnorm(nrow*ncol);
  Temp.attr("dim") = Dimension(nrow, ncol);
  
  NumericMatrix eib = as<NumericMatrix>(Temp);

  for(i = 0; i < nrow; ++i)
  { mean_eib(i) = mean(eib(i,_));
    sd_eib(i) = sd(eib(i,_));
    out(i)=0;
    for (j = 0; j < ncol; ++j){
      Tbij(i,j) = (eib(i,j)-mean_eib(i))/sd_eib(i);
      if (NumericMatrix::is_na(W(i,j))==FALSE){
        out(i)+=W(i,j)/mi(i)+ci(i)*Tbij(i,j)*W(i,j);
      }
    }

  }
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
