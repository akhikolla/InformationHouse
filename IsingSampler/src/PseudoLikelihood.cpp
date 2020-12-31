#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double PseudoLikelihood(NumericMatrix x, NumericMatrix graph, NumericVector thresholds, double beta, IntegerVector responses, bool logis){
  int n = x.nrow();
  int k = x.ncol();
  
  double logPS = 0;
  double e = 0;
  
  // for every response, compute log PL and sum:
  for (int i=0;i<n;i++){
    // for every node, compute log PL and sum:
    for (int j=0;j<k;j++){
      e=thresholds[j];
      for (int jj=0;jj<k;jj++){
        e += x(i,jj) * graph(j, jj);
      }
      logPS += x(i,j) * e  -  log(exp(responses[0]*e) + exp(responses[1]*e));
    }
  }
  
  if (!logis){
    return(exp(logPS));
  } else {
     return(logPS); 
  }
}
