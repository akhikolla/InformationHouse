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

#include "AsymMixedPdfIntegrandEvaluator.h"

typedef AsymMixedPdfIntegrandEvaluator AMPIE;

AMPIE::AsymMixedPdfIntegrandEvaluator(arma::vec eigP): eigenP(eigP) {}

std::complex<double> AMPIE::integrand(double x, double t, double maxError) {
  if (t == 0) {
    return x / (2 * M_PI);
  }
  std::complex<double> val;
  std::complex<double> I(0, 1);

  std::complex<double> sum = 0;
  std::complex<double> v(0, 12.0 * (-2.0 * t) / (M_PI * M_PI));
  double precision = std::pow(static_cast<double>(10), -15);
  for(int i = 0; i < eigenP.size(); i++) {
    if (std::fabs(eigenP[i]) > precision) {
      int sign = getSinhSign((v * eigenP[i]).imag());
      std::complex<double> sinhProdVal = sinhProd(v * eigenP[i], 1);
      if (sinhProdVal.imag() * sign <= 0) {
        sinhProdVal *= -1;
      }
      sum += std::log(sinhProdVal);
    }
  }
  return 1 / (2 * M_PI) * std::exp(sum) * std::exp(-I * t * x);
}
