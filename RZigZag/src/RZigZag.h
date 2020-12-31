// RZigZag.h : implements Zig-Zag and other PDMP samplers
//
// Copyright (C) 2017--2019 Joris Bierkens
//
// This file is part of RZigZag.
//
// RZigZag is free software: you can redistribute it and/or modify it
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
// along with RZigZag.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __RZIGZAG_H
#define __RZIGZAG_H

#define __INCLUDE_EIGEN
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

#include "ZigZag.h"
#include "LogisticSampler.h"
#include "GaussianSampler.h"
#include "VariousSamplers.h"
#include "PostProcess.h"

List SkeletonToList(const Skeleton& skel);

Skeleton ListToSkeleton(const List& listSkeleton);

#endif
