/*
 Estimation Volume Under ROC Surface.
 Author: Khanh To Duc
 */

#include <Rcpp.h>
using namespace Rcpp;

static inline double indvus(double a, double b, double c) {
  if((a < b) && (b < c)) {
    return 1.0;
  }else if((a < b) && (b == c)){
    return 0.5;
  }else if((a == b) && (b < c)){
    return 0.5;
  }else if((a == b) && (b == c)){
    return 1.0/6;
  }
  else{
    return 0.0;
  }
}

// [[Rcpp::export]]
double vusC(NumericVector tt, NumericMatrix dd){
  int nn = tt.size();
  double I_ijk = 0.0;
  double den = 0.0, num = 0.0, temp = 0.0;
  for(int i = 0; i < nn; i++){
    for(int j = 0; j < nn; j++){
      if(j != i){
        for(int k = 0; k < nn; k++){
          if((k != j) && (k != i)){
            I_ijk = indvus(tt[i], tt[j], tt[k]);
            temp = dd(i, 0)*dd(j, 1)*dd(k, 2);
            num += I_ijk*temp;
            den += temp;
          }
        }
      }
    }
  }
  double out = num/den;
  return out;
}

// [[Rcpp::export]]
double vusC_full(NumericVector tt1, NumericVector tt2, NumericVector tt3){
  int nn1 = tt1.size(), nn2 = tt2.size(), nn3 = tt3.size();
  double num = 0.0;
  for(int i = 0; i < nn1; i++){
    for(int j = 0; j < nn2; j++){
      for(int k = 0; k < nn3; k++){
        num += indvus(tt1[i], tt2[j], tt3[k]);
      }
    }
  }
  double out = num/(nn1*nn2*nn3);
  return out;
}

// [[Rcpp::export]]
List vusC_full_core(NumericVector tt1, NumericVector tt2, NumericVector tt3){
  int nn1 = tt1.size(), nn2 = tt2.size(), nn3 = tt3.size();
  NumericVector out1(nn1), out2(nn2), out3(nn3);
  NumericMatrix I_ij(nn2, nn1);
  NumericMatrix I_ik(nn3, nn1);
  for(int i = 0; i < nn1; i++){
    NumericMatrix M(nn3, nn2);
    for(int j = 0; j < nn2; j++){
      for(int k = 0; k < nn3; k++){
        M(k,j) = indvus(tt1[i], tt2[j], tt3[k]);
      }
    }
    NumericVector temp = colSums(M);
    I_ij(_, i) = temp;
    I_ik(_, i) = rowSums(M);
    out1[i] = sum(temp);
  }
  out2 = rowSums(I_ij);
  out3 = rowSums(I_ik);
  List out;
  out["ind1"] = out1;
  out["ind2"] = out2;
  out["ind3"] = out3;
  return out;
}

