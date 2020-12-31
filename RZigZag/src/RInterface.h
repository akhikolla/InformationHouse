// RInterface.h
//
// Copyright (C) 2017--2019 Joris Bierkens
//
// This file is part of RZigZag.
//
// This is an attempt to bring all zigzag functionality into a pure C++ setting
// with the idea of possible interfacing to both Python and R
//
// ZigZag is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RZigZag is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ZigZag.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __RINTERFACE_H
#define __RINTERFACE_H

#define __INCLUDE_EIGEN
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

VectorXd getUniforms(const long n);
VectorXd getStandardNormals(const long n);

static std::ostream& messageStream = Rcout;

#endif
