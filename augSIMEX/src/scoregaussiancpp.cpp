#include <Rcpp.h>
using namespace Rcpp;



NumericMatrix scoregaussian(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericVector weight, NumericVector offset) {
  int i,j;
  int nrow=DataM.nrow();
  int ncol=DataM.ncol();
  NumericMatrix out(nrow,ncol+1);
  double eta;
  
  for(i = 0; i < nrow; ++i)
  { eta=0;
    for (j = 0; j < ncol; ++j){
      eta += beta(j)*DataM(i,j);
    }
    eta+=offset(i);
    for (j = 0; j < ncol; ++j){
      out(i,j) = weight(i)*(Y(i)-eta)*DataM(i,j)/pow(beta(ncol),2);
    }
    out(i,ncol) = weight(i) * (- 1/beta(ncol) + pow(Y(i)-eta,2)/pow(beta(ncol),3));
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector scoregaussiancpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, NumericVector offset) {
  int i,j;
  int nrow=DataM.nrow();
  int ncol=DataM.ncol();
  NumericMatrix scorez0(nrow,ncol+1),scorez1(nrow,ncol+1);
  NumericVector value(ncol+1);
  
  scorez0=scoregaussian(beta,Y,DataM0,weight,offset);
  scorez1=scoregaussian(beta,Y,DataM1,weight,offset);
  
  for (j = 0; j < ncol+1; ++j){
    value(j)=0;
    for(i = 0; i < nrow; ++i){
      value(j)+=((1-DataM(i,ncol-1))*(scorez0(i,j)*(1-phat(i))-scorez1(i,j)*qhat(i))-DataM(i,ncol-1)*(scorez0(i,j)*phat(i)-scorez1(i,j)*(1-qhat(i))))/(1-phat(i)-qhat(i));
    }
  }
  
  return value;
}


/*** R

*/
