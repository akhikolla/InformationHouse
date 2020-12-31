// LogisticSampler.cpp
//
// Copyright (C) 2017--2019 Joris Bierkens
//
// This file is part of RZigZag.
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

#include "LogisticSampler.h"

double LogisticData::potential(const VectorXd& position) const {
  double val = 0;
  for (int j = 0; j < n_observations; ++j) {
    double innerproduct = position.dot(dataX.row(j));
    val += log(1 + exp(innerproduct)) - dataY(j) * innerproduct;
  }
  return val;
}

VectorXd LogisticData::gradient(const VectorXd& position) const {
  VectorXd grad(VectorXd::Zero(dim));
  for (int j = 0; j < n_observations; ++j) {
    double val = exp(dataX.row(j).dot(position));
    grad += dataX.row(j).transpose() * (val/(1+val) - dataY(j));
  }
  return grad;
}

MatrixXd LogisticData::hessian(const VectorXd& position) const {
  MatrixXd hess(MatrixXd::Zero(dim,dim));
  for (int j = 0; j < n_observations; ++j) {
    double innerproduct = position.dot(dataX.row(j));
    hess += (dataX.row(j).transpose() * dataX.row(j))* exp(innerproduct)/((1+exp(innerproduct)*(1+exp(innerproduct))));
  }
  return hess;
}

double LogisticData::getDerivative(const VectorXd& position, const int index) const {
  
  double derivative = 0;
  for (int j = 0; j < n_observations; ++j) {
    double val = exp(dataX.row(j).dot(position));
    derivative += dataX(j,index) * (val/(1+val) - dataY(j));
  }
  return derivative;
}

double LogisticData::getSubsampledDerivative(const VectorXd& position, const int index, const VectorXd& x_ref) const {
  int J = floor(n_observations*getUniforms(1)(0)); // randomly select observation
  if (x_ref.size() == 0)
    return n_observations * dataX(J,index) * (1.0/(1.0+exp(-dataX.row(J).dot(position))) - dataY(J));
  else
     return n_observations * dataX(J,index) * (1.0/(1.0+exp(-dataX.row(J).dot(position))) - 1.0/(1.0+exp(-dataX.row(J).dot(x_ref))));
}

MatrixXd LogisticData::dominatingHessian() const {
  
  MatrixXd domHessian(MatrixXd::Zero(dim,dim));
  for (int j = 0; j < n_observations; ++j) {
    domHessian += 0.25 * (dataX.row(j).transpose() * dataX.row(j));
  }
  return domHessian;
}

VectorXd LogisticData::getUniformBound() const {
  const VectorXd norms (dataX.rowwise().norm());
  VectorXd bounds(dim);
  
  for (int k = 0; k < dim; ++k) {
    bounds(k) = 0.0;
    for (int l = 0; l < n_observations; ++l) {
      double val = fabs(dataX(l,k) * norms(l));
      if (bounds(k) < val)
        bounds(k) = val;
    }
  }
  return 0.25 * bounds * n_observations;
}

void LogisticZZ::Initialize() {
  a = state.v.array() * data.gradient(state.x).array();
  b = (long double) sqrt((long double) dim) * data.dominatingHessian().rowwise().norm();
}
  
double LogisticZZ::getTrueIntensity() {
  return data.getDerivative(state.x, proposedEvent) * state.v(proposedEvent);
}

void LogisticZZ::updateBound() {
  if (acceptProposedEvent)
    a(proposedEvent) = -trueIntensity;
  else
    a(proposedEvent) = trueIntensity;
}

void LogisticCVZZ::Initialize() {
  if (x_ref.size() == 0) {
    x_ref = VectorXd::Zero(dim);
    grad_ref = newton(x_ref, data);
  }
  else {
    grad_ref = data.gradient(x_ref);
  }
  C = data.getUniformBound();
  b = (long double)sqrt((long double)dim) * C;
  updateBound();
}

double LogisticCVZZ::getTrueIntensity() {
  return state.v(proposedEvent)*(grad_ref(proposedEvent) + data.getSubsampledDerivative(state.x, proposedEvent, x_ref));
}

void LogisticCVZZ::updateBound() {
  a_ref = (state.v.cwiseProduct(grad_ref)).unaryExpr(&pospart);
  a = (state.x-x_ref).norm() * C + a_ref;
}
