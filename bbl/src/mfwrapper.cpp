#include <Rcpp.h>
#include <gsl/gsl_vector.h>
#include "bfgs.h"
using namespace Rcpp;

// [[Rcpp::export]]
List mfwrapper(NumericMatrix xi, IntegerVector weights, NumericMatrix qJ,
               IntegerVector Lv, NumericVector Eps, NumericVector priorCount){
  
  int n = xi.nrow();
  int m = xi.ncol();
  std::vector<short> L;
  std::vector<std::vector<short> > sv(n);
  for(int i=0; i<m; i++){
    L.push_back(Lv[i]);
    for(int k=0; k<n; k++)
      sv[k].push_back(xi(k,i));
  }
  
  std::vector<int> frq(n);
  for(int k=0; k<n; k++)
    frq[k] = weights(k);
  
  int mL=L.size();
  std::vector<std::vector<double> > hp(mL);
  std::vector<std::vector<std::vector<double> > > Jp(mL);
  
  double eps = Eps[0];
  double lkl = 0;
  double lnz = 0;
  double priorcount = priorCount(0);
  invC(sv, frq, L, lkl, lnz, hp, Jp, eps, priorcount);
  
  std::vector<std::vector<double> > h(m);
  std::vector<std::vector<std::vector<double> > > J(m);
  
  int ic=0;
  for(int i=0; i<m; i++){
    J[i].resize(m);
    h[i]=hp[ic];
    int jc=0;
    for(int j=0; j<m; j++){
      J[i][j] = Jp[ic][jc++];
      if(!qJ(i,j))
        std::fill(J[i][j].begin(), J[i][j].end(), 0.0); // non-interacting pairs
    }
    ic++;
  }
    
  List x = List::create(Named("h") = h, Named("J") = J, Named("lkh") = lkl, 
                        Named("lz")= lnz);
  return x;
}
