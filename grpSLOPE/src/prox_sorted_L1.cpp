// SLOPE: Sorted L1 Penalized Estimation (SLOPE)
// Copyright (C) 2015 Malgorzata Bogdan, Ewout van den Berg, Chiara Sabatti,
// Weijie Su, Emmanuel Candes, Evan Patterson
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
#include "proxSortedL1.h"
}

// [[Rcpp::export]]
NumericVector prox_sorted_L1_C(NumericVector y, NumericVector lambda) {
  size_t n = y.size();
  NumericVector x(n);
  evaluateProx(y.begin(), lambda.begin(), x.begin(), n, NULL);
  return x;
}
