#include <Rcpp.h>
using namespace Rcpp;


// function for rolling mean
// [[Rcpp::export]]
NumericVector rollMeanC(NumericVector x, int k, char type) {

  int n = x.size();
  int n_k = n-k+1;

  /*
  // only used für leading or trailing missings (not yet implemented)
  int i_start,i_end;
  if(type=='c'){
    i_start = ceil((double) k/2)-1;
    i_end = n - i_start + 1;
  }else if(type=='r'){
    i_start = 0;
    i_end = n-k+1;
  }else if(type=='l'){
    i_start = k-1;
    i_end = n;
  }

  return(i_start);
  */
  NumericVector x_out(n_k);


  for(int i=0;i<n_k;i++){
    x_out[i] = mean(x[seq(i,i+k-1)]);
  }

  return x_out;
}


// function for rolling sum
// [[Rcpp::export]]
NumericVector rollSumC(NumericVector x, int k, char type) {

  int n = x.size();
  int n_k = n-k+1;

  /*
  // only used für leading or trailing missings (not yet implemented)
  int i_start,i_end;
  if(type=='c'){
  i_start = ceil((double) k/2)-1;
  i_end = n - i_start + 1;
  }else if(type=='r'){
  i_start = 0;
  i_end = n-k+1;
  }else if(type=='l'){
  i_start = k-1;
  i_end = n;
  }

  return(i_start);
  */
  NumericVector x_out(n_k);


  for(int i=0;i<n_k;i++){
    x_out[i] = sum(x[seq(i,i+k-1)]);
  }

  return x_out;
}


// function for weighted ratio
//' @name PointEstimates
//' @title Weighted Point Estimates
//'
//' @description Predefined functions for weighted point estimates in package `surveysd`.
//'
//' @param x numeric vector
//' @param w weight vector
//'
//' @details Predefined functions are weighted ratio and weighted sum.
//'
//' @return
//' Each of the functions return a single numeric value
//' @examples
//' x <- 1:10
//' w <- 10:1
//' weightedRatio(x,w)
//' @export
// [[Rcpp::export]]
double weightedRatio(NumericVector x, NumericVector w) {

  int n = x.size();
  double upper=0;
  double lower=0;

  for(int i=0;i<n;i++){
    if(x[i]==1){
      upper = upper+w[i];
    }
    if(!NumericVector::is_na(x[i])){
      lower = lower+w[i];
    }
  }
  double out=upper/lower*100;

  return out;
}

// function for weighted sum
//' @rdname PointEstimates
//' @examples
//' x <- 1:10
//' w <- 10:1
//' weightedSum(x,w)
//' @export
// [[Rcpp::export]]
double weightedSum(NumericVector x, NumericVector w) {

  int n = x.size();
  double out=0;

  for(int i=0;i<n;i++){
    if(!NumericVector::is_na(x[i])){
      out = out+x[i]*w[i];
    }
  }

  return out;
}

