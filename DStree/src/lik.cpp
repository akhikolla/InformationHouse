#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

 double lik(NumericVector S, NumericVector pi,
             NumericVector ncens, NumericVector nuncens) {
    int n = S.size();
    double out=0;
    for(int i=0;i<n;i++){
    if(S[i]==0)
        S[i]=1;
	if(pi[i]==0)
        pi[i]=1;
	 
	 out+= nuncens[i]*log(pi[i])+ncens[i]*log(S[i]);
	 }
	 return(out);
	 }
	 
