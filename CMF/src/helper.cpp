#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
void p_gradUsparse(NumericMatrix Xm, NumericMatrix Gm, NumericMatrix CUm, NumericMatrix OUm, NumericMatrix Cm, NumericVector I, NumericVector T, NumericVector Rowm, NumericVector Colm){
 int N = Xm.nrow(), K = Gm.ncol();
 int row, r, c, k;
 double temp;
 if(I(0)==1) {
   for(row=0;row<N;row++) {
     r = Xm(row,0) - 1; c = Xm(row,1) - 1;
     temp = 0;
     for(k=0;k<K;k++) {
       temp = temp + CUm(r,k)*OUm(c,k);
     }
     temp = temp - Xm(row,2) + Rowm(r) + Colm(c);
     for(k=0;k<K;k++) {
       Gm(r,k) = Gm(r,k) + T(0) * (temp * OUm(c,k) + CUm(r,k)*Cm(c,k));
     }
   }
 } else {
   for(row=0;row<N;row++) {
     r = Xm(row,0) - 1; c = Xm(row,1) - 1;
     temp = 0;
     for(k=0;k<K;k++) {
       temp = temp + CUm(c,k)*OUm(r,k);
     }
     temp = temp - Xm(row,2) + Rowm(r) + Colm(c);
     for(k=0;k<K;k++) {
       Gm(c,k) = Gm(c,k) + T(0) * (temp * OUm(r,k) + CUm(c,k)*Cm(r,k));
     }
   }
 }
 //return(Gm);
}


// [[Rcpp::export]]
void p_updatePseudoData(NumericMatrix Xm, NumericMatrix U1m, NumericMatrix U2m, NumericVector Rv, NumericVector Cv){
  int N = Xm.nrow(), K = U1m.ncol();
  int row, r, c, k;
  double temp;
  for(row=0;row<N;row++) {
    r = Xm(row,0) - 1; c = Xm(row,1) - 1;
    temp = 0;
    for(k=0;k<K;k++) {
      temp = temp + U1m(r,k)*U2m(c,k);
    }
    temp = temp + Rv(r) + Cv(c);
    Xm(row,2) = temp;
  }
  //return(Xm);
}

// [[Rcpp::export]]
NumericVector p_updateTau(NumericMatrix Xm, NumericMatrix U1m, NumericMatrix U2m, NumericMatrix cov1m, NumericMatrix cov2m, NumericVector Rv, NumericVector Cv, NumericVector nu1v, NumericVector nu2v){
  int N = Xm.nrow(), K = U1m.ncol();
  int row, r, c, k;
  double temp;
  NumericVector out(1);

  for(row=0;row<N;row++) {
    r = Xm(row,0) - 1; c = Xm(row,1) - 1;
    temp = 0;
    for(k=0;k<K;k++) {
      temp = temp + U1m(r,k)*U2m(c,k);
    }
    temp = temp + Rv(r) + Cv(c);
    temp = Xm(row,2) - temp;
    temp = temp*temp;
    for(k=0;k<K;k++) {
      temp = temp + cov1m(r,k)*U2m(c,k)*U2m(c,k) + U1m(r,k)*U1m(r,k)*cov2m(c,k) + cov1m(r,k)*cov2m(c,k);
    }
    temp = temp + nu1v(r) + nu2v(c);
    out(0) = out(0) + temp;
  }
  return(out);
}

// [[Rcpp::export]]
List p_updateMean(NumericMatrix Xm, NumericMatrix U1m, NumericMatrix U2m, NumericVector I, NumericVector Mv){
  int N = Xm.nrow(), K = U1m.ncol();
  int row, r, c, k;
  double temp;
  NumericVector Nv, Cv;

  if(I(0)==1) {
    Nv = NumericVector(U1m.nrow());
    Cv = NumericVector(U1m.nrow());
  } else {
    Nv = NumericVector(U2m.nrow());
    Cv = NumericVector(U2m.nrow());
  }

  for(row=0;row<N;row++) {
    r = Xm(row,0) - 1; c = Xm(row,1) - 1;
    temp = 0;
    for(k=0;k<K;k++) {

      temp = temp + U1m(r,k)*U2m(c,k);
    }
    if(I(0)==1) {
      temp = Xm(row,2) - temp - Mv(c);
      Nv(r) = Nv(r) + temp;
      Cv(r)++;

    } else {
      temp = Xm(row,2) - temp - Mv(r);
      Nv(c) = Nv(c) + temp;
      Cv(c)++;
    }
  }
  List ret;
  ret["sum"] = Nv;
  ret["count"] = Cv;
  return(ret);
}


// [[Rcpp::export]]
void p_covUsparse(
  NumericMatrix Xm, NumericMatrix Cm, NumericMatrix OUm, NumericMatrix OCm, NumericVector I, NumericVector T){
  int N = Xm.nrow(), K = Cm.ncol();
  int row, r, c, k;

  if(I(0)==1) {
    for(row=0;row<N;row++) {
      r = Xm(row,0) - 1; c = Xm(row,1) - 1;
      for(k=0;k<K;k++) {
        Cm(r,k) = Cm(r,k) + T(0) * (OUm(c,k)*OUm(c,k) + OCm(c,k));
      }
     }
  } else {
    for(row=0;row<N;row++) {
      r = Xm(row,0) - 1; c = Xm(row,1) - 1;
      for(k=0;k<K;k++) {
        Cm(c,k) = Cm(c,k) + T(0)* (OUm(r,k)*OUm(r,k) + OCm(r,k));
      }
    }
  }
  //return(Cm);
}



