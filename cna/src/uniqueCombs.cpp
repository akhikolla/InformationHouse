#include <Rcpp.h>
using namespace Rcpp;

// Count the number of unique values in each column
// [[Rcpp::export]]
IntegerVector C_countUniques(const IntegerMatrix x){
  int n=x.nrow(), p=x.ncol();
  IntegerVector nVals(p);
  for (int j=0; j<p; j++){
    int cnt=n - sum(as<IntegerVector>(duplicated(x(_, j))));
//    Rcpp::Rcout << "[Rcpp::Rcout] " << cnt << std::endl;
    nVals[j]=cnt;
    }
  return nVals;
}


// Identify duplicated rows in a matrix
// [[Rcpp::export]]
LogicalVector C_duplicatedMat(const IntegerMatrix x){
  int n=x.nrow(), p=x.ncol();
  IntegerVector nv=C_countUniques(x);
  IntegerVector f(nv.size());
  for (int i=0; i<p; i++){
    if (i==0){ 
      f(i)=1;
      //Rcpp::Rcout << "[Rcpp::Rcout] " << f(i) << std::endl;
    } else {
      f(i)=f(i-1)*nv(i-1);
      //Rcpp::Rcout << "[Rcpp::Rcout] " << f(i) << std::endl;
    }
  }  
  IntegerVector rowConfig(n);
  for (int i=0; i<n; i++){
    for (int j=0; j<p; j++){
      rowConfig(i) += (x(i, j)-1) * f(j);
    }  
  }
  return duplicated(rowConfig);
}
  
// Eliminate duplicated rows in a matrix
// [[Rcpp::export]]
IntegerMatrix C_uniqueMat(const IntegerMatrix x){
  LogicalVector keep=!C_duplicatedMat(x);
  int n=x.nrow(), p=x.ncol(), m=sum(as<IntegerVector>(keep));
  IntegerMatrix out(m, p);
  int ii=0;
  for (int i=0; i<n; i++){
    if (keep(i)){
      for (int j=0; j<p; j++){
        out(ii, j)=x(i,j);
      }    
    ii++;  
    }  
  }
  return out;
}  
    
// Subsetting an IntegerMatrix by column
// [[Rcpp::export]]
IntegerMatrix C_selectCols(const IntegerMatrix x, const IntegerVector idx){
  int n=x.nrow(), r=idx.size();
  IntegerMatrix out(n, r);
  for (int i=0; i<n; i++){
    for (int j=0; j<r; j++){
      out(i, j)=x(i, idx[j]-1);
    }
  }
  return out;
}

// Unique rows from a subset of the columns of a matrix
// [[Rcpp::export]]
IntegerMatrix C_uniqueCombs(const IntegerMatrix x, const IntegerVector idx){
  return C_uniqueMat(C_selectCols(x, idx));
}


/* R

expand <- function(x){
  stopifnot(x>0, x %% 1 == 0)
  C_expand(x)
}

xx <- expand(2:4)

C_duplicatedMat(xx)
C_uniqueMat(xx)

set.seed(12)
xs <- xx[sample(nrow(xx), 12, replace = TRUE), ]
identical(unique(xs), C_uniqueMat(xs))
identical(xx[, c(1, 3)], C_selectCols(xx, c(1, 3)))

idx <- c(1, 3)
identical(unique(xs[, idx]), C_uniqueCombs(xs, idx))

library(microbenchmark)
microbenchmark(unique(xs[, idx]), C_uniqueCombs(xs, idx))
*/

  
