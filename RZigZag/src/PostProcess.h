// PostProcess.h
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

#ifndef __POSTPROCESS_H
#define __POSTPROCESS_H

#include "ZigZag.h"


class PostProcess {
public:
  
  PostProcess(Skeleton& skel): skel{skel}, covarianceEstimated{false}, asvarEstimated{false}, meansEstimated{false} {};

  // postprocessing functions
  void estimateCovariance(const SizeType coordinate = -1, bool zeroMeans = false);
  void estimateAsymptoticVariance(const int n_batches = 100, const SizeType coordinate = -1, bool zeroMeans = false);
  void estimateESS(const int n_batches = 100, const SizeType coordinate = -1, bool zeroMeans = false);
  void estimateMoment(const int p, const SizeType coordinate = -1);
  void sample(const int n_samples, const SizeType coordinate = -1);

  VectorXd getESS() const { return ess;}
  VectorXd getAsVar() const { return asVar;}
  VectorXd getMoment() const { return moment;}
  VectorXd getMeans() const { return means;}
  MatrixXd getCovarianceMatrix() const { return covarianceMatrix; }
  MatrixXd getSamples() const { return samples;}
  VectorXd getSampleTimes() const { return sampleTimes;}
  
private:
  const Skeleton& skel;
  bool covarianceEstimated;
  bool asvarEstimated;
  bool meansEstimated;
  
  MatrixXd covarianceMatrix;
  VectorXd means;
  VectorXd asVar;
  VectorXd ess;
  VectorXd moment;
  MatrixXd samples;
  VectorXd sampleTimes;

};

#endif
