// PostProcess.cpp
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

#include "PostProcess.h"


void PostProcess::estimateAsymptoticVariance(const int n_batches, const SizeType coordinate, bool zeroMeans) {
  if (n_batches <= 0)
    throw std::range_error("n_batches should be positive.");
  const SizeType dim = (coordinate < 0 ? skel.Positions.rows() : 1);
  const double t_max = skel.Times[skel.size-1];
  const double batch_length = t_max / n_batches;
  
  double t0 = skel.Times[0];
  VectorXd x0(dim), x1(dim);
  if (coordinate < 0)
    x0 = skel.Positions.col(0);
  else
    x0 = skel.Positions.row(coordinate).col(0);
  
  MatrixXd batchMeans(dim, n_batches);
  
  int batchNr = 0;
  double t_intermediate = batch_length;
  VectorXd currentBatchMean = VectorXd::Zero(dim);
  
  for (int i = 1; i < skel.size; ++i) {
    double t1 = skel.Times[i];
    if (coordinate < 0)
      x1 = skel.Positions.col(i);
    else
      x1 = skel.Positions.row(coordinate).col(i);
    
    while (batchNr < n_batches - 1 && t1 > t_intermediate) {
      VectorXd x_intermediate = x0 + (t_intermediate - t0) / (t1 - t0) * (x1 - x0);
      batchMeans.col(batchNr) = currentBatchMean + (t_intermediate - t0) * (x_intermediate + x0)/(2 * batch_length);
      
      // initialize next batch
      currentBatchMean = VectorXd::Zero(dim);
      batchNr++;
      t0 = t_intermediate;
      x0 = x_intermediate;
      t_intermediate = batch_length * (batchNr + 1);
    }
    currentBatchMean += (t1 - t0) * (x1 + x0)/(2 * batch_length);
    t0 = t1;
    x0 = x1;
  }
  batchMeans.col(batchNr) = currentBatchMean;
  
  if (!zeroMeans)
    means = batchMeans.rowwise().sum()/n_batches;
  else
    means = VectorXd::Zero(dim);
  
  MatrixXd meanZeroBatchMeans = batchMeans.colwise() - means;
  asVar = batch_length * meanZeroBatchMeans.rowwise().squaredNorm()/(n_batches - 1);
  
  asvarEstimated = true;
  // ESS = (covarianceMatrix.diagonal().array()/asVarEst.array() * t_max).matrix();
}

void PostProcess::estimateESS(const int n_batches, const SizeType coordinate, bool zeroMeans) {
  
  if (!covarianceEstimated)
    estimateCovariance(coordinate, zeroMeans);
  if (!asvarEstimated)
    estimateAsymptoticVariance(n_batches, coordinate, zeroMeans);
  double t_max = skel.Times[skel.size-1];
  ess = covarianceMatrix.diagonal().array()/asVar.array() * t_max;
}

void PostProcess::estimateMoment(const int p, const SizeType coordinate) {
  
  const SizeType dim = (coordinate < 0 ? skel.Positions.rows() : 1);
  
  const double t_max = skel.Times[skel.size-1];
  
  double t0 = skel.Times[0];
  VectorXd x0(dim), x1(dim);
  if (coordinate < 0)
    x0 = skel.Positions.col(0);
  else
    x0 = skel.Positions.row(coordinate).col(0);
  
  ArrayXd momentArray = ArrayXd::Zero(dim);

  for (int i = 1; i < skel.size; ++i) {
    double t1 = skel.Times[i];
    if (coordinate < 0)
      x1 = skel.Positions.col(i);
    else
      x1 = skel.Positions.row(coordinate).col(i);
    // the following expression equals \int_{t_0}^{t_1} x(t) (x(t))^T d t
    momentArray += (t1 - t0) * (x1.array().pow(p+1) - x0.array().pow(p+1))/((p+1)*t_max * (x1-x0).array());
    t0 = t1;
    x0 = x1;
  }
  moment = VectorXd(momentArray);
}

void PostProcess::estimateCovariance(const SizeType coordinate, bool zeroMeans) {
  
  const SizeType dim = ( coordinate < 0 ? skel.Positions.rows() : 1);
  
  const double t_max = skel.Times[skel.size-1];
  
  double t0 = skel.Times[0];
  VectorXd x0(dim), x1(dim);
  if (coordinate < 0)
    x0 = skel.Positions.col(0);
  else
    x0 = skel.Positions.row(coordinate).col(0);
  
  covarianceMatrix = MatrixXd::Zero(dim, dim);
  means = VectorXd::Zero(dim);
  
  for (int i = 1; i < skel.size; ++i) {
    double t1 = skel.Times[i];
    if (coordinate < 0)
      x1 = skel.Positions.col(i);
    else
      x1 = skel.Positions.row(coordinate).col(i);
    // the following expression equals \int_{t_0}^{t_1} x(t) (x(t))^T d t
    covarianceMatrix += (t1 - t0) * (2 * x0 * x0.transpose() + x0 * x1.transpose() + x1 * x0.transpose() + 2 * x1 * x1.transpose())/(6 * t_max);
    if (!zeroMeans)
      means += (t1 - t0) * (x1 + x0) /(2 * t_max);
    t0 = t1;
    x0 = x1;
  }
  covarianceMatrix -= means * means.transpose();
  covarianceEstimated = true;

}

void PostProcess::sample(const int n_samples, const SizeType coordinate) {
  
  const SizeType dim = ( coordinate < 0 ? skel.Positions.rows() : 1);
  const int n_steps = skel.Times.size();
  if (n_steps < 2)
    throw "Skeleton::sample: skeleton size < 2.";
  const double t_max = skel.Times(n_steps-1);
  const double dt = t_max / (n_samples+1);
  
  double t_current = dt;
  double t0 = skel.Times(0);
  double t1;
  VectorXd x0(dim), x1(dim);
  if (coordinate < 0)
    x0 = skel.Positions.col(0);
  else
    x0 = skel.Positions.row(coordinate).col(0);

  samples = MatrixXd(dim, n_samples);
  sampleTimes = VectorXd(n_samples);
  int n_sampled = 0; // number of samples collected
  
  for (int i = 1; i < n_steps; ++i) {
    if (coordinate < 0)
      x1 = skel.Positions.col(i);
    else
      x1 = skel.Positions.row(coordinate).col(i);
    t1 = skel.Times(i);
    while (t_current < t1 && n_sampled < n_samples) {
      samples.col(n_sampled) = x0 + (x1-x0) * (t_current - t0)/(t1-t0);
      ++n_sampled;
      t_current = t_current + dt;
    }
    x0 = x1;
    t0 = t1;
  }
  
  for (int i=0; i < n_samples; ++i) {
    sampleTimes[i] = (i+1)*dt;
  }

}
