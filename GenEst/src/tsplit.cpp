#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calcRateC(NumericMatrix M, NumericMatrix Aj,
  NumericVector days, NumericMatrix searches){
// M = (estimated) mortalities (ncarc x nsim)
// Aj = (simulated) arrival intervals (ncarc x nsim)
// days = vector of all search days (on all units) = c(0, ...)
// searches = (ncarc x ndays) array of 0s & 1s; indicates search days for carcs
  const int nsim = M.ncol();
  const int x = M.nrow();
  const int nsearch = days.size() - 1;
  int xi, i, i0, i1;
// di = index into days for carcass xi [so c(1, 0, 1, 1, 0) -> c(0, 2, 3, 0, 0)]
  int **di;
  di = new int *[x];
  for (xi = 0; xi < x; xi++) {
    di[xi] = new int [nsearch + 1]();
  }
  for (xi = 0; xi < x; xi++){
    i = 0;
    for (int si = 0; si < days.size(); si++){
      if (int(searches(xi, si)) == 1){
        di[xi][i] = si;
        i++;
      }
    }
  }
  NumericMatrix rate(nsim, nsearch);
  double tmprate;
  for (int simi = 0; simi < nsim; simi++){
    for (int xi = 0; xi < x; xi++){
      i0 = int(Aj(xi, simi)) - 1;
      i1 = i0 + 1;
      tmprate = M(xi, simi)/(days[di[xi][i1]] - days[di[xi][i0]]);
      for (i = di[xi][i0]; i < di[xi][i1]; i++){
        rate(simi, i) += tmprate; 
      }
    }
  }
// need to prorate the rate for the unmonitored parts of the season [not done]
  return(rate);
}

// [[Rcpp::export]]
NumericMatrix calcTsplitC(NumericMatrix rate, NumericVector days, NumericVector times){
  // rate = array of rates
  // days = master search schedule with all search days
  // times = vector of bin boundaries for temporal split
  const int nsim = rate.nrow();
  const int ntimes = times.size();
  NumericMatrix splits(ntimes - 1, nsim);
  int ti, si;
  double *Ft;
  Ft = new double[ntimes]();
  for (int simi = 0; simi < nsim; simi++){
    for (ti = 1; ti < ntimes; ti++){
      Ft[ti] = 0;
      si = 1;
      while (si < days.size()){
        if(days[si] <= times[ti]){
          Ft[ti] = Ft[ti] + (days[si] - days[si - 1]) * rate(simi, si - 1);
        } else {
          Ft[ti] = Ft[ti] + (times[ti] - days[si - 1]) * rate(simi, si - 1);
          break;
        }
        si++;
      }
      splits(ti - 1, simi) = (Ft[ti] - Ft[ti - 1]);
    }
  }
  return(splits);
}
