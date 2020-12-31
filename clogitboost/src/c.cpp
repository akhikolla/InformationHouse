#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

double score(NumericVector g, NumericVector y, NumericVector w, int n, NumericVector group, double theta){
	int i;
	int j;
  double temp1 = 0.0, temp2 = 1.0;
  int samegroup = 0;
	double output = 0.0;
	for(i = 0; i < n; i++) {
	  temp1 = 0;
	  temp2 = 0;
	  if (y[i] == 1){
      for(j = 0; j < n; j++) {
	      if (group[i] == group[j]){
		      samegroup = 1;
		    } else {
		      samegroup = 0;
		    }
	      temp1 += g[j] * exp(w[j] + (theta) * g[j]) * samegroup;
	      temp2 += exp(w[j] + (theta) * g[j]) * samegroup;
      }
      output += g[i] * y[i];
      if (temp2 != 0 ){
        output -= temp1 / temp2;
      }
	  }
  }
	return output;
}



// [[Rcpp::export]]
double seconddev(NumericVector g, NumericVector y, NumericVector w, int n, NumericVector group, double theta) {
  int i;
	int j;
  double temp1 = 0.0, temp2 = 0.0, temp1d = 0.0;
  int samegroup = 0;
  double output = 0.0;
  
  for(i = 0; i < n; i++) {
	  temp1 = 0;
	  temp2 = 0;
	  temp1d = 0;
    if (y[i] == 1){
	    for(j = 0; j < n; j++) {
        if (group[i] == group[j]){
          samegroup = 1;
        }	else {
		      samegroup = 0;
		    }
        temp1 += g[j] * exp(w[j] + (theta) * g[j]) * samegroup;
	      temp2 += exp(w[j] + (theta) * g[j]) * samegroup;
	      temp1d += g[j] * g[j] * exp(w[j] + (theta) * g[j]) * samegroup;
      }
	
	    if (temp2 !=0 ){
      output -= (temp1d * temp2 - temp1 * temp1) / temp2 / temp2;
		  }
	  }
  }
    
	return output;
}


// [[Rcpp::export]]
NumericVector persamplegrad(NumericVector fx, NumericVector y, int n, NumericVector group) {
  int i;
	int j;
  double temp2 = 1.0;
	int samegroup = 0;
	NumericVector output(n);
	
  for(i = 0; i < n; i++) {
	  temp2 = 0;  
	  for(j = 0; j < n; j++) {
	    if (group[i] == group[j]){
		    samegroup = 1;
		  }	else {
		    samegroup=0;
		  }
      temp2 += exp(fx[j]) * samegroup;	
    }
		output[i] = y[i] - exp(fx[i]) / temp2;
	}
	return output;
}

// [[Rcpp::export]]
double likelihood(NumericVector fx, NumericVector y, int n, NumericVector group) {
  int i;
	int j;
  double temp2 = 1.0;
	int samegroup = 0;
	double output=0.0;

  for(i = 0; i < n; i++) {
		if (y[i] == 1){
	    temp2 = 0;	    
	    for(j = 0; j < n; j++) {	
	      if (group[i] == group[j]){
		      samegroup = 1;
		    }	else {
		      samegroup = 0;
        }	
	      temp2 += exp(fx[j]) * samegroup;	
      }
	    output += fx[i] * y[i] - log(temp2);
	  }
	}
  return output;
}
	
	
