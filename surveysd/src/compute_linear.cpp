#include <Rcpp.h>
using namespace Rcpp;

//' @rdname computeFrac
//' @export
// [[Rcpp::export]]
NumericVector computeLinear(double curValue,
                            double target,
                            const NumericVector& x,
                            const NumericVector& w,
                            double boundLinear = 10) {
  NumericVector f(x.size());
  if(x.size()!=w.size()){
    stop("x and w of different length!");
  }
  if(x.size()==1){
    f[0] = w[0]/curValue*target;
    return f;
  }else{
    double h = 0.0;
    double j = 0.0;
    double N = 0.0;

    for(int i = 0; i < x.size(); i++){
      h += w[i]*x[i];
      j += w[i]*x[i]*x[i];
      N += w[i];
    }

    double b = (target-N*j/h)/(h-N*j/h);
    double a = (N-b*N)/h;


    for(int i = 0; i < x.size(); i++)
      f[i] = a*x[i] + b;

    //apply bounds
    for(int i = 0; i < x.size(); i++){
      if (f[i] < 1.0/boundLinear)
        f[i] = 1.0/boundLinear;
      if (f[i] > boundLinear)
        f[i] = boundLinear;
    }

    return f;
  }
}

//' @rdname computeFrac
//' @export
// [[Rcpp::export]]
NumericVector computeLinearG1(double curValue,
                              double target,
                              const NumericVector& x,
                              const NumericVector& w,
                              double boundLinear = 10) {
  NumericVector f(x.size());
  f=computeLinear(curValue,target,x,w,boundLinear);
  for(int i = 0; i < x.size(); i++){
    if (f[i]*w[i] < 1.0){
      f[i]=1/w[i];
    }
  }
  return f;
}
