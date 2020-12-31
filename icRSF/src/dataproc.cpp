#include <Rcpp.h>
using namespace Rcpp;

// Functions for processing data for likelihood functions

// [[Rcpp::export]]
NumericMatrix dmat(NumericVector id, NumericVector time, IntegerVector result, double phi1, double phi0, double negpred) {
  NumericVector utime = unique(time);
  utime.sort();
  int J = utime.size(), nsub = unique(id).size(), nobs = id.size(), i, j;
  NumericVector rid(nobs);
  NumericMatrix Cm(nsub, J + 1), Dm(nsub, J + 1);
  utime.push_back(utime[J-1] + 1); 
  //initiate with 1 for Cm
  std::fill(Cm.begin(), Cm.end(), 1);
  //Calculate Cm  
  j = 0;
  for (i = 1; i < nobs; i++) {
    if (id[i] != id[i-1]) j++;
    rid[i] = j;
  }
  for (j = 0; j <= J; j++) {
    for (i = 0; i < nobs; i++) {
      if (result[i]==0) {
        if (time[i] >= utime[j]) {
          Cm(rid[i], j) *= 1 - phi1; 
        } else {
          Cm(rid[i], j) *= phi0;
        }
      } else {
        if (time[i] >= utime[j]) {
          Cm(rid[i], j) *= phi1; 
        } else {
          Cm(rid[i], j) *= 1 - phi0;
        }      
      }
    }
  }
  //Calculate Dm
  for (i = 0; i < nsub; i++) {
    Dm(i, 0) = Cm(i, 0);
    for (j = 1; j < J+1; j++) {
    Dm(i, j) = negpred*(Cm(i, j) - Cm(i, j-1));
    }
  }
  return Dm; 
}

// [[Rcpp::export]]
IntegerVector getrids(NumericVector id, int nsub){
  IntegerVector rid(nsub);
  int j = 0;
  rid[0] = 1;
  for (int i = 1; i < id.size(); i++) {
    if (id[i] != id[i-1]) {
      j++;
      rid[j] = i + 1;
    }
  }
  return rid;
}

// Fucntion to convert time dependent coavriate matrix
IntegerVector timeIDX(NumericVector time, NumericVector utime) {
  int nobs = time.size(), J = utime.size();
  IntegerVector result(nobs);
  result.fill(-1);
  for (int i = 0; i < nobs; i++) {
    for (int j = 0; j < J; j++) {
      if (time[i] == utime[j]) {
        result[i] = j;
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix timeMat(int nsub, int J, NumericVector time, NumericVector utime, NumericMatrix Xmat) {
  IntegerVector timeid = timeIDX(time, utime);
  int next, nobs = timeid.size(), nbeta = Xmat.ncol(), i = 0;
  NumericMatrix result(nsub*J, nbeta);
  for (int j = 0; j < nobs; j++) {
    if (j == nobs - 1) {
      next = J;
    } else if (timeid[j+1] == 0) {
      next = J;
    } else {
      next = timeid[j + 1];
    }
    for (int k = timeid[j]; k < next; k++) {
      for (int m = 0; m < nbeta; m++) {
        result(i, m) = Xmat(j, m);
      }
      i++;
    }
  }
  return result;
}
