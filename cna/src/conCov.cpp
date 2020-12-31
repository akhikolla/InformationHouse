
#include <Rcpp.h>
using namespace Rcpp;


// Calculate con and cov of conjunctions
// [[Rcpp::export]]
NumericVector C_conj_conCov(const IntegerVector cols, const NumericMatrix x, 
                            const NumericVector y, const IntegerVector f){
  int n=x.nrow(), p=cols.size();
  NumericVector Sums(3), conCov(2);
  for (int i=0; i<n; i++){
    double minX=1;
    for (int j=0; j<p; j++){
      double Xval=x(i, cols[j]-1);
      if (Xval<minX) minX=Xval;
    };
    double Yval=y[i];
    double minXY=std::min(minX,Yval);
    Sums(0)+=minX*f(i);
    Sums(1)+=Yval*f(i);
    Sums(2)+=minXY*f(i);
  };
  conCov(0) = Sums(2)/Sums(0);
  conCov(1) = Sums(2)/Sums(1);
  return conCov;
}

// Calculate con and cov of disjunctions
// [[Rcpp::export]]
NumericVector C_disj_conCov(const IntegerVector cols, const NumericMatrix x, 
                            const NumericVector y, const IntegerVector f){
  int n=x.nrow(), p=cols.size();
  NumericVector Sums(3), conCov(2);
  for (int i=0; i<n; i++){
    double maxX=0;
    for (int j=0; j<p; j++){
      double Xval=x(i, cols[j]-1);
      if (Xval>maxX) maxX=Xval;
    };
    double Yval=y[i];
    double XY=std::min(maxX,Yval);
    Sums(0)+=maxX*f(i);
    Sums(1)+=Yval*f(i);
    Sums(2)+=XY*f(i);
  };
  conCov(0) = Sums(2)/Sums(0);
  conCov(1) = Sums(2)/Sums(1);
  return conCov;
}

/* R
set.seed(1234)
x <- matrix(runif(20), 4)
cols <- c(2, 4)
y <- runif(nrow(x))
f1 <- rep(1L, length(y))
f <- seq_along(x)
x
y
f
extractColumns(x, cols)
M(x, cols)
M2(x, cols, y)
C_conj_conCov(x, cols, y, f)

*/
