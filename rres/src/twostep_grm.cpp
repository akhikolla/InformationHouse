#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix twostep_grm_cpp(NumericMatrix genotype, NumericVector freq, NumericMatrix est0) {
  int nind = genotype.nrow();

  NumericMatrix est(nind, nind);
  
  for(int i = 0; i < nind; i++){
    for(int j = i; j < nind; j++){
      NumericVector al = 2.0 * freq / (1.0 + est0(i, j)) + est0(i, j) / (1.0 + est0(i, j));
      NumericVector Vl_inv = 1.0 / ((pow(al - 2.0 * freq, 2) + est0(i,j) * pow(al - 1.0, 2) 
                                       - est0(i,j)/2.0) / 4.0 / freq / (1.0 - freq) + est0(i,j) * (1.0 - est0(i,j)/2.0)/2.0 + 0.25);
      NumericVector wl = Vl_inv / sum(Vl_inv);
      
      est(i, j) = sum(wl * (genotype(i,_) * genotype(j,_) - al*(genotype(i,_) + genotype(j,_)) + 
        4.0 * freq * (al - freq)) / freq / (1.0 - freq)) / 2.0;
      if(i != j){
        est(j, i) = est(i, j);
      }
    }
  }
  
  return est;
}

