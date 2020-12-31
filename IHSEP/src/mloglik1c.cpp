#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double mloglik1c(NumericVector jtms, double TT, Function nu,
		 NumericVector gcoef, Function Inu){
  double res=0.0;
  double s=0.0;
  int nj=jtms.size();
  NumericVector lmbd(nj),exeff(nj);
  lmbd=nu(jtms);
  exeff[0]=0.0;
  for(int j=1;j<nj;j++){
    exeff[j]=exeff[j-1]*exp(-gcoef[1]*(jtms[j]-jtms[j-1]))+
      gcoef[0]*exp(-gcoef[1]*(jtms[j]-jtms[j-1]));
  }
  for(int j=0;j<nj;j++)lmbd[j]+=exeff[j];
  for(int j=0;j<nj;j++) res -= log(lmbd[j]);
  // Rcpp::Rcout<<res<<std::endl;
  res += as<double>(Inu(TT));
  // Rcpp::Rcout<<res<<std::endl;
  for(int j=0;j<nj;j++)s += 1-exp(-gcoef[1]*(TT-jtms[j]));
  s *= gcoef[0]/gcoef[1];
  // Rcpp::Rcout<<"s="<<s<<std::endl;
  res += s;
  return(res);
}
