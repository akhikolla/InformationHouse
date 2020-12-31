#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

//' This function summarizes the z matrix into 2 K x S matrices
//' where K is the number of communities and S is the number of species.
//' The matrix "res1" holds the values for dat=1
//' The matrix "res0" holds the values for dat=0
// [[Rcpp::export]]
Rcpp::List getks(IntegerMatrix z, int ncommun, IntegerMatrix dat) {
  int nspp=z.ncol();
  int nlinhas=z.nrow();
  IntegerMatrix res1(ncommun,nspp);
  IntegerMatrix res0(ncommun,nspp);

  for(int i=0; i<nlinhas;i++){
    for (int j=0; j<nspp; j++){
      if (dat(i,j)==1) {
        res1(z(i,j)-1,j)=res1(z(i,j)-1,j)+1;
      }
      if (dat(i,j)==0){
        res0(z(i,j)-1,j)=res0(z(i,j)-1,j)+1;
      }
    }
  }

  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("nks1") = res1,
                                          Rcpp::Named("nks0") = res0);
  return(resTemp);
}

//' This function summarizes the z matrix into a L x K matrix
//' where K is the number of communities and L is the number of locations
// [[Rcpp::export]]
IntegerMatrix getlk(IntegerMatrix z, IntegerVector locid, int ncommun, int nloc) {
  int nlinhas=z.nrow();
  int nspp=z.ncol();
  IntegerMatrix res(nloc,ncommun);

  for(int i=0; i<nlinhas;i++){
    for (int j=0; j<nspp; j++){
      res(locid[i]-1,z(i,j)-1)=res(locid[i]-1,z(i,j)-1)+1;
    }
  }

  return(res);
}

// This function helps with multinomial draws
int whichLessDVPresenceFast(double value, NumericVector prob) {
  int res=-1;
  double probcum = 0;

  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples z's
// [[Rcpp::export]]
IntegerMatrix samplez(NumericMatrix ltheta, NumericMatrix l1minustheta, NumericMatrix lphi,
                      IntegerMatrix dat1, IntegerVector locid,
                      NumericMatrix randu, int ncommun, int nloc) {

  IntegerMatrix zmat(dat1.nrow(),dat1.ncol());
  NumericVector prob(ncommun);
  int znew;

  for(int i=0; i<dat1.nrow();i++){
    for (int j=0; j<dat1.ncol(); j++){
      for (int k=0; k<ncommun; k++){
        prob(k)=ltheta(locid(i)-1,k)+dat1(i,j)*lphi(k,j)+(1-dat1(i,j))*l1minustheta(k,j);
      }
      prob=prob-max(prob);
      prob=exp(prob);
      prob=prob/sum(prob);

      //multinomial draw
      znew=whichLessDVPresenceFast(randu(i,j),prob);
      zmat(i,j)=znew+1;
    }
  }
  return zmat;
}

//' This function converts vmat into theta
// [[Rcpp::export]]
NumericMatrix convertVtoTheta(NumericMatrix vmat,
                              NumericVector prod) {
  NumericMatrix res(vmat.nrow(),vmat.ncol());

  for(int j=0; j<vmat.ncol();j++){
    res(_,j)=vmat(_,j)*prod;
    prod=prod*(1-vmat(_,j));
  }

  return (res);
}
