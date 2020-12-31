#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector latentX_calc_cpp(NumericVector parAlpha, NumericVector epsilonT, int ndays) {
  
  NumericVector latentX(ndays);
  
  latentX[0] = epsilonT[0]; // latentX on the first day
  // AR(1) calculation for latentX
  for(int i_day = 1; i_day < ndays; ++i_day) {
    latentX[i_day] = (parAlpha[i_day] * latentX[i_day-1]) + epsilonT[i_day];
  }
  
  return latentX;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# test R code
*/
