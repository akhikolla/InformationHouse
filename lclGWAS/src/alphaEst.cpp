/*****************************/
/* Jiaxing Lin               */
/* 05-13-2016                :w*/
/* Estimate Alpha without    */
/*****************************/

#include "Rcpp.h"

using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
List alphaEst(
    Rcpp::NumericVector dtime,
    Rcpp::NumericVector delta) 
{
  double d = 0;
  double r = 0;
  
  int m    = max(dtime); // number of time breaks for alphas

  Rcpp::NumericVector alpha(m);  
  int NPat = dtime.size();// number of patient	
  for(int j=0; j < m; j++){
  	d = 0;
    r = 0;
    for(int i=0; i < NPat; i++)
  	{
    	if(dtime[i] == j+1 && delta[i] == 1)
      		d++;
    	if(dtime[i] > j+1)
      		r++; 	
  	} 	
  	alpha[j] = r/(r+d);
  }
  return List::create(
    Named("alphaEst")
    = alpha
  );
  
}




