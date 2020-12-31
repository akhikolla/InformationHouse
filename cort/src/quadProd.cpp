#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double quadProd(const NumericMatrix a,
                const NumericMatrix b,
                const NumericVector kern,
                const NumericMatrix other_a,
                const NumericMatrix other_b,
                const NumericVector other_kern) {

  int d = a.ncol();
  int n = a.nrow();
  int other_n = other_a.nrow();
  double temp;
  double rez = 0.0;

  for (int i=0; i<n; i++){
    for (int j=0; j<other_n; j++){
      temp = kern(i)*other_kern(j);
      if(temp != 0.0){
        for (int dim=0; dim<d; dim++){
          temp *= std::max(std::min(b(i,dim),other_b(j,dim)) - std::max(a(i,dim),other_a(j,dim)),0.0);
          if(temp == 0.0){
            break;
          }
        }
        rez += temp;
      }
    }
  }
  return(rez);
}
