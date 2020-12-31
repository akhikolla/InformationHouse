#include<Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector GoFsim(int nsim, int n, Function fun) {
  NumericVector res(nsim);
  NumericVector V(n);
  for(int i=0;i<nsim;i++){
    GetRNGstate();
    V =-log(runif(n));
    PutRNGstate();
    res[i]=as<double>(fun(V));

  }
 // NumericVector temp=wrap(res);
  std::sort(res.begin(),res.end());
  return res;
  }
