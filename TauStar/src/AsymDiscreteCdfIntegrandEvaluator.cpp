/***
 * Copyright (C) 2016 Luca Weihs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "AsymDiscreteCdfIntegrandEvaluator.h"

typedef AsymDiscreteCdfIntegrandEvaluator ADCIE;

ADCIE::AsymDiscreteCdfIntegrandEvaluator(
  arma::vec eigP, arma::vec eigQ): eigenP(eigP), eigenQ(eigQ) {}

std::complex<double> ADCIE::integrand(double x, double t, double maxError) {
  if (t == 0) {
    return x / (2 * M_PI);
  }
  std::complex<double> val;
  std::complex<double> I(0, 1);

  std::complex<double> sum = 0;
  for(int i = 0; i < eigenP.size(); i++) {
    for(int j = 0; j < eigenQ.size(); j++) {
      sum += -0.5 * std::log(1.0 - 8.0 * I * t * eigenP[i] * eigenQ[j]);
    }
  }
  return 1 / (2 * M_PI) * std::exp(sum) * (1.0 - std::exp(-I * t * x)) / (I * t);
}
