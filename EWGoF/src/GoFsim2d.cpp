#include<Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix GoFsim2d(int nsim, int n, Function fun1, Function fun2) {
  NumericMatrix res(2,nsim);
  NumericVector V(n);
  for(int i=0;i<nsim;i++){
    GetRNGstate();
    V = -log(runif(n));
    PutRNGstate();
    res(1,i)=as<double>(fun1(V));
    res(2,i)=as<double>(fun2(V));
  }
//NumericMatrix temp=wrap(res);
  std::sort(res.begin(),res.end());
  return res;
  }
