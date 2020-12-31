#include <Rcpp.h>
using namespace Rcpp;

// rankEnsCpp compute ranking of ensembles of forecasts

//    Copyright (C) 2016 MeteoSwiss

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//   This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.



// [[Rcpp::export]]
NumericVector rankEnsCpp(NumericMatrix ens) {
  int n = ens.nrow(), nens = ens.ncol() ;
  NumericVector ranks(n, 1.0) ; 
  for (int i = 1; i < n; i++){
    for (int j = 0; j < i; j++){
      NumericVector enstmp1 = ens(i, _);
      NumericVector enstmp2 = ens(j, _);
      NumericVector enstmp(2*nens), ensranks(2*nens);
      std::copy(enstmp1.begin(), enstmp1.end(), enstmp.begin());
      std::copy(enstmp2.begin(), enstmp2.end(), enstmp.begin() + nens);
      ensranks = match(enstmp, clone(enstmp).sort()) ;
      int s1 = 0, s2 = 0;
      for (int k = 0; k < nens; k++){
        s1 += ensranks(k);
        s2 += ensranks(k + nens);
      }
      if (s1 > s2){
        ranks(i) += 1.0;
      } else if (s1 < s2) {
        ranks(j) += 1.0;
      } else {
        ranks(i) += 0.5;
        ranks(j) += 0.5;
      }
    }
  }
  return(ranks);
}

