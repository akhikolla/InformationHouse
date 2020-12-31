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

#ifndef TauStar_AsymDiscretePdfIntegrandEvaluator
#define TauStar_AsymDiscretePdfIntegrandEvaluator

#include "IntegrandEvaluator.h"
#include "RcppArmadillo.h"

class AsymDiscretePdfIntegrandEvaluator : public IntegrandEvaluator {
protected:
  arma::vec eigenP;
  arma::vec eigenQ;
public:
  AsymDiscretePdfIntegrandEvaluator(arma::vec eigP, arma::vec eigQ);
  std::complex<double> integrand(double x, double t, double maxError);
};

#endif
