#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double mloglik1d(NumericVector jtms, double TT, Function nu,
		 NumericVector gcoef, Function Inu){
  int d=gcoef.size()/2;
  int nj=jtms.size();
  NumericVector lmbd(nj);
  NumericVector exeff(d,0.0);
  lmbd=nu(jtms);
  double res=0, s=0.0;
  int j,k;
  res -= log(lmbd[0]);
  j=1;
  while(j<nj){
    for(k=0;k<d;k++){
      exeff[k] += gcoef[k];
      exeff[k] *= exp(-gcoef[d+k]*(jtms[j]-jtms[j-1]));
      lmbd[j] += exeff[k];
    }
    res -= log(lmbd[j]);
    j++;
  }
  // Rcpp::Rcout<<res<<std::endl;
  res += as<double>(Inu(TT));
  // Rcpp::Rcout<<res<<std::endl;
  for(k=0;k<d;k++){
    s=0.0;
    for(j=0;j<nj;j++){
      s += (1-exp(-gcoef[k+d]*(TT-jtms[j])))*gcoef[k]/gcoef[k+d];
    }
    // s*=gcoef[k]/gcoef[k+d];
    // Rcpp::Rcout<<"s="<<s<<std::endl;
    res += s;
  }
  return(res);
}
