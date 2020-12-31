#include <Rcpp.h>
#include <gsl/gsl_vector.h>
#include "bfgs.h"
using namespace Rcpp;

// [[Rcpp::export]]
List pseudo_mle(NumericMatrix xi, IntegerVector freq,
                LogicalMatrix qJ, IntegerVector Lv, NumericVector Lambda, 
                IntegerVector Nprint, IntegerVector Itmax, NumericVector Tol,
                LogicalVector Naive, IntegerVector Verbose, 
                LogicalVector Lzhalf){
  
  int n = xi.nrow();
  int m = xi.ncol();
  std::vector<short> L(m);
  std::vector<std::vector<short> > sv(n);
  for(int k=0; k<n; k++) for(int i=0; i<m; i++)
    sv[k].push_back(xi(k,i));

  std::vector<int> frq(n);
  for(int k=0; k<n; k++)
    frq[k] = freq(k);
  
  int nbad=0;
  std::vector<bool> bad(m);
  std::vector<std::vector<bool> > qj(m);
  
  for(int i=0; i<m; i++){
    for(int j=0; j<m; j++)
      qj[i].push_back(qJ(i,j));
    short xmin=xi(0,i);
    short xmax=xi(0,i);
    for(int k=1; k<n; k++){ 
      if(xmax<xi(k,i)) xmax=xi(k,i);
      if(xmin>xi(k,i)) xmin=xi(k,i);
    }
    if(xmax==xmin){ 
      bad[i]=true;
      nbad++;
    }
    else bad[i]=false;
    L[i]=Lv(i);
  }
  if(nbad>0){
    std::vector<short> v(m);
    for(int i=0; i<m; i++){      // add an extra instance to regularize
      if(!bad[i]) v[i]=0;
      else if(L[i]>0) v[i]=0;
      else{
        v[i]=1;
        L[i]=1;
      }
    }
    sv.push_back(v);
    frq.push_back(1);
    n++;
  }
  
  std::vector<std::vector<double> > h(m);
  std::vector<std::vector<std::vector<double> > > J(m);

  double lambda = Lambda[0];
  double lambdah = Lambda[1];
  int nprint = Nprint[0];
  unsigned int Imax = Itmax[0];
  double tol = Tol[0];
  int verbose = Verbose[0];
  bool naive = Naive[0];
  bool lzhalf = Lzhalf[0];
  
  double lks = 0;
  double lz = 0;
  for(int i0=0; i0<m; i0++){
    double z=0;
    bool failed=false;
    lks += lpr_psl(i0, sv, frq, qj[i0], L, lambda, lambdah, 
                   h[i0], J[i0], nprint, Imax, tol, verbose, z, naive, 
                   failed, lzhalf);
    if(failed)
      Rcpp::Rcerr << " Warning: failed to converge in pseudo\n";
    lz += z;
    Rcpp::checkUserInterrupt();
  }
  lks /= n;
  
  List x = List::create(Named("h") = h, Named("J") = J, 
                        Named("lkh") = lks, Named("lz") = lz, Named("L") = L);
  return x;
}

// [[Rcpp::export]]
NumericVector predict_class(IntegerVector xid, IntegerVector Ly, List h, List J,
                            NumericVector lz, NumericVector py, 
                            LogicalVector Naive){

  int m = xid.length();
  int ly = Ly[0];
  
  NumericVector E(ly);
  for(int iy=0; iy<ly; iy++){
    double e=0;
    List hy = h[iy];
    for(int i=0; i<m; i++){
      if(xid(i)==0) continue;
      NumericVector hi = hy[i];
      if(hi.length() < xid(i)) continue;
      e += hi(xid(i)-1);
      if(Naive[0]) continue;
      List Jy = J[iy];
      List Ji = Jy[i];
      for(int j=0; j<m; j++){
        if(j==i || xid(j)==0) continue;
        NumericMatrix Jj = Ji[j];
        if(Jj.nrow()<xid(i) ||  Jj.ncol()<xid(j)) continue;
        e += Jj(xid(i)-1,xid(j)-1)/2.0;
      }
    }
    E[iy] = e - lz[iy] + log(py[iy]);
    Rcpp::checkUserInterrupt();
  }
  
  return E;
}
