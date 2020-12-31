#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//' comp_surv_cpp

// [[Rcpp::export]]
int comp_surv_cpp(double time1, double event1, double time0, double event0) {
  int comp = 0;
  if (((event1==1) & (time0>time1)) | ((event1==1) & (event0==0) & (time1==time0))){
    comp = -1;
  } else if (((event0==1) & (time1>time0)) | ((event1==0) & (event0==1) & (time1==time0))){
    comp = 1;
  }
  return comp;
}

//' mat_comp_surv_cpp

// [[Rcpp::export]]
NumericMatrix mat_comp_surv_cpp(NumericVector time1, NumericVector event1, NumericVector time0, NumericVector event0){
  int n1 = time1.length();
  int n0 = time0.length();
  NumericMatrix M(n1, n0);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n0; j++) {
      M(i, j) = comp_surv_cpp(time1[i], event1[i], time0[j], event0[j]);
    }
  }
  return M;
}

//' comp_cont_cpp

// [[Rcpp::export]]
int comp_cont_cpp(double event1, double event0, std::string direction) {
  int comp = 0;
  if (direction == ">"){
    if (event0 < event1){
      comp = -1;
    } else if (event0 > event1){
      comp = 1;
    }
  } else if (direction == "<"){
    if (event0 > event1){
      comp = -1;
    } else if (event0 < event1){
      comp = 1 ;
    }
  }
  return comp;
}

//' mat_comp_cont_cpp

// [[Rcpp::export]]
NumericMatrix mat_comp_cont_cpp(NumericVector event1, NumericVector event0, std::string direction){
  int n1 = event1.length();
  int n0 = event0.length();
  NumericMatrix M(n1, n0);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n0; j++) {
      M(i, j) = comp_cont_cpp(event1[i], event0[j], direction);
    }
  }
  return M;
}

//' sign_cpp

// [[Rcpp::export]]
int sign_cpp(int x){
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}

//' comp_repeated_cpp

// [[Rcpp::export]]
int comp_repeated_cpp(NumericVector time1, NumericVector event1, double fu1, NumericVector time0, NumericVector event0, double fu0){
  int comp = 0;
  int n = time1.length();
  if (fu1 < fu0){
    for (int k = 0; k < n; k++) {
      if (time0[k] > fu1){
        time0[k] = fu1;
        event0[k] = 0;
      }
    }
  } else if (fu1 > fu0){
    for (int k = 0; k < n; k++) {
      if (time1[k] > fu0){
        time1[k] = fu0;
        event1[k] = 0;
      }
    }
  }
  if (sum(event1) > sum(event0)) {
    comp = -1;
  } else if (sum(event1) < sum(event0)) {
    comp = 1;
  } else {
    for (int k = 0; k < n; k++) {
      comp = comp + comp_surv_cpp(time1[k], event1[k], time0[k], event0[k]); 
    }
    comp = sign_cpp(comp);
  }
  return comp;
}

//' mat_comp_repeated_cpp

// [[Rcpp::export]]
NumericMatrix mat_comp_repeated_cpp(NumericMatrix time1, NumericMatrix event1, NumericVector fu1, NumericMatrix time0, NumericMatrix event0, NumericVector fu0){
  int n1 = fu1.length();
  int n0 = fu0.length();
  NumericMatrix M(n1, n0);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n0; j++) {
      NumericVector t1 = time1(i,_);
      NumericVector e1 = event1(i,_);
      NumericVector t0 = time0(j,_);
      NumericVector e0 = event0(j,_);
      M(i, j) = comp_repeated_cpp(t1, e1, fu1[i], t0, e0, fu0[j]);
    }
  }
  return M;
}
