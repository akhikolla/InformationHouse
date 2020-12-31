#include <Rcpp.h>
using namespace Rcpp;

// EnsRoca.R compute area under the ROC curve

//    Copyright (C) 2016 MeteoSwiss

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

// [[Rcpp::export]]
NumericVector rankCpp(NumericVector x) {
  int n=x.size();
  NumericVector ranks(n);
  ranks = match(x, clone(x).sort());
  NumericVector uniqueranks = unique(ranks);
  int nunique=uniqueranks.size();
  if (nunique != n) {
    for (int k = 0; k < nunique; k++){
      double nind = sum(ranks == uniqueranks(k));
      std::replace(ranks.begin(), ranks.end(), uniqueranks(k), uniqueranks(k) + (nind - 1.0)/2.0);
    }
  }
  return(ranks); 
}

// [[Rcpp::export]]
NumericVector EnsRocaCpp(NumericMatrix ens, NumericMatrix obs) {
  int n = ens.nrow(), ncat = ens.ncol();
  double ntotal = ens.nrow();
  NumericVector rocarea(ncat);
  for (int i = 0; i < n; i++){
    double probsum = sum(ens(i, _));
    if (probsum != 1.0) {
      ens(i, _) = ens(i, _) / probsum;
    }
  }
  for (int j = 0; j < ncat; j++){
    NumericVector enscat = ens( _, j);
    NumericVector obscat = obs(_, j);
    double nevent = sum(obscat);
    NumericVector ensrank(ntotal);
    ensrank = rankCpp(enscat);
    double meanrank = sum(ensrank * obscat) / nevent;
    rocarea[j] = (meanrank - (nevent + 1.0)/2.0) / (ntotal - nevent);
    if (nevent == 0) {
      rocarea[j] = NA_REAL ;
    }
  }
  return(rocarea) ; 
} 
