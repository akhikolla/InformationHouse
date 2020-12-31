#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix scorelogit(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericVector weight, NumericVector offset) {
  int i,j;
  int nrow=DataM.nrow();
  int ncol=DataM.ncol();
  NumericMatrix out(nrow,ncol);
  double eta=0;
  double expeta;
  
  for(i = 0; i < nrow; ++i)
  { eta=0;
    for (j = 0; j < ncol; ++j){
      eta += beta(j)*DataM(i,j);
    }
    eta+=offset(i);
    expeta=exp(eta);
    for (j = 0; j < ncol; ++j){
      out(i,j) = weight(i)*Y(i)*DataM(i,j)-weight(i)*expeta/(1+expeta)*DataM(i,j);
    }
  }
  
  return out;
}


// [[Rcpp::export]]
NumericVector scorelogitcpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, NumericVector offset) {
  int i,j;
  int nrow=DataM.nrow();
  int ncol=DataM.ncol();
  NumericMatrix scorez0(nrow,ncol),scorez1(nrow,ncol);
  NumericVector value(ncol);
  
  scorez0=scorelogit(beta,Y,DataM0,weight,offset);
  scorez1=scorelogit(beta,Y,DataM1,weight,offset);
  
  for (j = 0; j < ncol; ++j){
    value(j)=0;
    for(i = 0; i < nrow; ++i){
      value(j)+=((1-DataM(i,ncol-1))*(scorez0(i,j)*(1-phat(i))-scorez1(i,j)*qhat(i))-DataM(i,ncol-1)*(scorez0(i,j)*phat(i)-scorez1(i,j)*(1-qhat(i))))/(1-phat(i)-qhat(i));
    }
  }
  
  return value;
}


/*** R

*/
