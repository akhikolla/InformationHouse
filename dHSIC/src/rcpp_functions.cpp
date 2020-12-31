#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix gaussian_grammat_rcpp(NumericMatrix x, double bandwidth, int n, int d){
  double xnorm = 0;
  NumericMatrix K(n,n);
  for(int i = 0; i < n; ++i){
    int j = i;
    while(j < n){
      for(int l = 0; l < d; ++l){
        xnorm += pow(x(i,l)-x(j,l), 2.0);
      }
      K(i,j) = exp(-xnorm/(2.0*pow(bandwidth, 2.0)));
      K(j,i) = K(i,j);
      xnorm = 0.0;
      ++j;
    }
  }
  return K;
}


// [[Rcpp::export]]
NumericMatrix discrete_grammat_rcpp(NumericMatrix x, int n, int d){
  int tmp = 0;
  NumericMatrix K(n,n);
  for(int i = 0; i < n; ++i){
    int j = i;
    while(j < n){
      for(int l = 0; l < d; ++l){
        tmp += (x(i,l)==x(j,l));
      }
      K(i,j) = double(tmp==d);
      K(j,i) = K(i,j);
      tmp = 0;
      ++j;
    }
  }
  return K;
}


// [[Rcpp::export]]
NumericMatrix shuffle_grammat_rcpp(NumericMatrix K, NumericVector perm, int n){
  NumericMatrix Kperm(n,n);
  for(int i = 0; i < n; ++i){
    int j = i;
    while(j < n){
      Kperm(i,j) = K(perm[i], perm[j]);
      Kperm(j,i) = Kperm(i,j);
      ++j;
    }
  }
  return Kperm;
}


// [[Rcpp::export]]
double median_bandwidth_rcpp(NumericMatrix x, int n, int d){
  int len = n;
  if(n > 1000){
    len = 1000;
  }
  double xnorm = 0.0;
  int lentot = len*(len+1)/2-len;
  int middle = lentot/2;
  NumericVector bandvec(lentot);
  int count = 0;
  for(int i = 0; i < len; ++i){
    int j = i+1;
    while(j < len){
      for(int l = 0; l < d; ++l){
        xnorm += pow(x(i,l)-x(j,l), 2.0);
      }
      bandvec[count] = xnorm;
      xnorm = 0.0;
      ++j;
      ++count;
    }
  }
  NumericVector v = clone(bandvec);
  std::nth_element(v.begin(), v.begin() + middle, v.end());
  double bandwidth = v[middle];
  bandwidth = sqrt(bandwidth*0.5);
  return bandwidth;
}
