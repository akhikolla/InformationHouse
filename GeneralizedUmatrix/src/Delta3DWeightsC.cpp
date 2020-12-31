#include <RcppArmadillo.h>
//using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube Delta3DWeightsC(Rcpp::NumericVector vx,Rcpp::NumericVector Datasample) {
  
  Rcpp::IntegerVector x_dims = vx.attr("dim");
  arma::cube x(vx.begin(), x_dims[0], x_dims[1], x_dims[2], false);
  //std::cout<<Datasample.length()<<" "<<x.n_slices<<std::endl;
  //std::cout<<x.slice(0)(1,1)<<std::endl;
  //std::cout<<x(1,1,2)<<std::endl;
  //arma::mat result(x.n_rows, x.n_cols);;
  for (unsigned int i = 0; i < x.n_slices; i++) {
    //std::fill(result.begin(),result.end(),Datasample(i));
    //x.slice(i)=x.slice(i)-result;
    x.slice(i)=x.slice(i)-Datasample(i);
//result.col(i) = arma::conv_to<arma::colvec>::from(arma::mean(x.slice(i)));  
  }
  
  return x;
}
