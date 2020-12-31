# include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector d_gauss_cpp(NumericVector X, double x, double sigma){
  NumericVector dis(X.length());  
  
  for(int i = 0; i < X.length(); i++){
    dis(i) = 2. / sigma * (X(i) - x);
  }
  return(dis);
}

// [[Rcpp::export]]
NumericVector d_mat52_cpp(NumericVector X, double x, double sigma){
  NumericVector s(X.length());
  double tmp;

  for(int i = 0; i < X.length(); i++){
    tmp = (x - X(i))/sigma;
    if(tmp > 0){
      s(i) = ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1. + sqrt(5.) * tmp + 5./3. * tmp * tmp);
    }else{
      if(tmp == 0){
        s(i) = 0;
      }else{
        tmp = std::abs(tmp);
        s(i) = -((10./3. - 5.) * tmp - 5. * sqrt(5.)/3. * tmp * tmp) / ((1. + sqrt(5.) * tmp + 5./3. * tmp * tmp));
      }
    }
  }
  return(s/sigma);
}

// [[Rcpp::export]]
NumericVector d_mat32_cpp(NumericVector X, double x, double sigma){
  NumericVector s(X.length());
  double tmp;

  for(int i = 0; i < X.length(); i++){
    tmp = (x - X(i))/sigma;
    if(tmp > 0){
      s(i) = -3*tmp / (1 + sqrt(3.) * tmp);
    }else{
      if(tmp == 0){
        s(i) = 0;
      }else{
        tmp = -tmp;
        s(i) = 3*tmp / (1 + sqrt(3.) * tmp);
      }
    }
  }
  return(s/sigma);
}
