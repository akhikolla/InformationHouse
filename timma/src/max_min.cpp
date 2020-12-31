#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//using namespace arma;

// max and min function for 3d array in Rcpp

// [[Rcpp::export]]
List maxcpp(NumericVector tau, int x, int y ,int j) {
   arma::cube mycube(tau.begin(), x, y, j);
   arma::vec sub(j);   
   arma::mat max(x,y), max_idx(x,y);
   for (int row = 0;row < x;row++){
     for (int col = 0;col < y;col++){
       for (int drug = 0;drug < j;drug++){
         sub(drug) = mycube(row, col, drug);
       }
       // get the index for max and min value
       arma::uword maxidx;
       // get max and min by slice
       max(row, col) = sub.max(maxidx);
       max_idx(row, col) = maxidx+1;
     }
   }
   return Rcpp::List::create(
                              Rcpp::Named("max") = max,
                              Rcpp::Named("max_idx") = max_idx
                              ) ;
}

// [[Rcpp::export]]
List mincpp(NumericVector tau, int x, int y ,int j) {
   arma::cube mycube(tau.begin(), x, y, j);
   arma::vec sub(j);
   arma::mat min(x,y), min_idx(x,y);
   arma::mat max(x,y), max_idx(x,y);
   for (int row = 0;row < x;row++){
     for (int col = 0;col < y;col++){
       for (int drug = 0;drug < j;drug++){
         sub(drug) = mycube(row, col, drug);
       }
       // get the index for max and min value
       arma::uword minidx;
       // get max and min by slice
       min(row, col) = sub.min(minidx);
       min_idx(row, col) = minidx+1;
     }
   }
   return Rcpp::List::create(
                              Rcpp::Named("min") = min,
                              Rcpp::Named("min_idx") = min_idx
                              ) ;
}


// [[Rcpp::export]]
arma::mat sumcpp(NumericVector tau, int x, int y ,int j) {
  arma::cube mycube(tau.begin(), x, y, j);
  arma::vec sub(j);
  arma::mat sum_cube(x,y);
  for(int row = 0;row < x;row++){
    for(int col = 0;col < y;col++){
      int count = 0;
      for(int drug = 0;drug<j;drug++){
        if (R_IsNA(mycube(row, col, drug))){
          sub(drug) = 0;
        }else{
          sub(drug) = mycube(row, col, drug);
          count++;
        }        
      }
      sum_cube(row, col) = sum(sub)/count;
    }
  }
  return(sum_cube);
}


// [[Rcpp::export]]
List maxcpp1(NumericVector tau, int x, int y) {
   arma::mat mymat(tau.begin(), x, y);
   arma::rowvec sub = max(mymat); 
   arma::rowvec row(y);
   arma::colvec col(x);
   for(int j = 0; j < y; j++){
     col = mymat.col(j);
     arma::uword maxidx;

     col.max(maxidx);
     row(j) = maxidx+1;
   }
   return Rcpp::List::create(
                              Rcpp::Named("max") = sub,
                              Rcpp::Named("max_idx") = row
                              ) ;
}


// [[Rcpp::export]]
List mincpp1(NumericVector tau, int x, int y) {
   arma::mat mymat(tau.begin(), x, y);
   arma::rowvec sub = min(mymat); 
   arma::rowvec row(y);
   arma::colvec col(x);
   for(int j = 0; j < y; j++){
     col = mymat.col(j);
     arma::uword minidx;

     col.min(minidx);
     row(j) = minidx+1;
   }
   return Rcpp::List::create(
                              Rcpp::Named("min") = sub,
                              Rcpp::Named("min_idx") = row
                              ) ;
}

// [[Rcpp::export]]
arma::rowvec sumcpp1(NumericVector tau, int x, int y) {
  arma::mat mymat(tau.begin(), x, y);
  arma::colvec sub(x);
  arma::rowvec s(y);
  for(int col = 0; col < y; col++){
    sub = mymat.col(col);
    int count = 0;
    for(int row = 0; row < x; row++){
      if(R_IsNA(sub(row))){
        sub(row) = 0;
      }else{
        
        count++;
      }
    }
    s(col) = sum(sub)/count;
  }

  return(s);
}
