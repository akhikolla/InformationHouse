// Optimized scatter
//
//'Get the matrix of optimized scatters used in locating the change points.
//'
//'@param II_ A D x N matrix where D is the maximum no. of segments (Kmax+1) and N is the no. of windows
//'@param X_ An N x r dataframe where N is the no. of windows and r the no. of running statistics monitored
//'@param H_ A D x N matrix where D is the maximum no. of segments (Kmax+1) and N is the no. of windows
//'
//'@return \item{II}{A matrix of optimized scatters}
//'@return \item{H}{A matrix of candidate changes point locations}

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List getScatterMatrix(NumericMatrix II_, NumericMatrix X_, NumericMatrix H_) {
  NumericMatrix X(X_);
  int M = X.ncol();
  int N = X.nrow();
  NumericMatrix u(N, N);
  NumericMatrix K(N, N);
  NumericMatrix fullCumSum(N+1, N+1);
  NumericVector diagCumSum(N+1);

  for (int i = 0; i < N; i++) {
    for(int j = i; j < N; j++) {
      double total= 0.0;
      for (int k = 0; k < M; k++) {
        double diff= X(i,k) - X(j,k);
        total += diff * diff;
      }
      u(i, j) = u(j, i) = total;
    }
  }

  double medianTimesTwo = 2 * median(u);

  // Upper triangle only!
  for (int i = 0; i < N; i++) {
    for (int j = 0; j <= i; j++) {
      K(j, i) = -u(j, i) / medianTimesTwo;
    }
  }

  // Upper triangle only!
  double top, left, leftTop;
  fullCumSum(0, N) = 0.0;
  for (int i = 0; i < N; i++) {
    top = 0.0;
    leftTop = 0.0;
    for (int j = 0; j <= i; j++) {
      left = fullCumSum(j + 1, i);
      fullCumSum(j + 1, i + 1) = top = left + top + exp(K(j, i)) - leftTop;
      leftTop = left;
    }
    fullCumSum(i + 1, i + 1) = fullCumSum(i, i + 1) + fullCumSum(i, i + 1) + exp(K(i, i)) - fullCumSum(i, i); // for left element, take top element on diagonal
  }

  for (int i = 0; i < N + 1; i++) {
    diagCumSum(i) = fullCumSum(i, i);
  }

  for (int j = 0; j < N; j++) {
    const double fullCumSum_j_j = diagCumSum(j + 1);
    for (int i = 0; i < j; i++) {
      const double temp = j - i + 1.0;
      K(i, j) = temp - (fullCumSum_j_j + diagCumSum(i) - 2.0 * fullCumSum(i, j + 1)) / temp;
    }
  }

  NumericMatrix II(II_);
  IntegerMatrix H(H_);
  int L = H.nrow();

  for (int i = 0; i < N; i++) {
    II(0, i) = K(0, i);
  }

  for(int k=1; k < L; ++k){
    for(int i=k; i < N; ++i){
      for(int j=k-1; j < i; ++j){
        double tmp = II(k-1, j) + K(j+1, i);  // Possible optimized criterion==> this is tried for all possible j<i
        if(tmp < II(k, i)){                  // By always retaining the least, we will end up with the minimum
          II(k, i) = tmp;
          H(k, i) = j+1;//fix indexing differences between R and C++
        }
      }
    }
  }

  return List::create(II, H);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#library(microbenchmark)
#library(testthat)

#S=1000
#X1=cbind(rnorm(S,0,1),rnorm(S,0,1)) #phase1: Means=0
#X2=cbind(rnorm(S,1,1),rnorm(S,1,1)) #phase2: Means=1
#X=rbind(X1,X2)
#X=scale(X) #each variable should be scaled to have an overall variance of 1
#microbenchmark(Y <- getScatterMatrix(X))
*/
