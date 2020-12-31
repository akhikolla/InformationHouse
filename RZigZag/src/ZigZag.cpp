// ZigZag.cpp
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

#include "ZigZag.h"

#define ACCURACY_CHECK 1e-6

Skeleton::Skeleton(const int dim, int initialSize) {
  if (initialSize < 1)
    initialSize = DEFAULTSIZE;
  Positions = MatrixXd(dim, initialSize);
  Velocities = MatrixXd(dim,initialSize);
  dimension = dim;
  Times = VectorXd(initialSize);
  capacity = initialSize;
  size = 0;
}

void Skeleton::Resize(const int factor) {
  capacity *= factor;
  Times.conservativeResize(capacity);
  Positions.conservativeResize(dimension, capacity);
  Velocities.conservativeResize(dimension, capacity);
}

void Skeleton::Push(const State& state, const double finalTime) {
  if (size >= capacity)
    Resize();
  Velocities.col(size) = state.v;
  if (finalTime < 0 || state.t < finalTime) {
    Times[size] = state.t;
    Positions.col(size) = state.x;
  }
  else {
    Times[size] = finalTime;
    double previousTime = Times[size-1];
    VectorXd previousPosition = Positions.col(size-1);
    Positions.col(size) = previousPosition + (finalTime - previousTime) * (state.x - previousPosition) / (state.t - previousTime);
  }
  size++;
}

void Skeleton::ShrinkToCurrentSize() {
  Times.conservativeResize(size);
  Positions.conservativeResize(dimension, size);
  Velocities.conservativeResize(dimension, size);
  capacity = size;
}

bool RejectionSampler::simulationStep() {
  // returns true if a switch is accepted
  
  acceptProposedEvent = false;
  proposeEvent(); // this moves the full sampler state ahead in time
  trueIntensity = getTrueIntensity();
  
  if (!rejectionFree) {
    double bound = getBound();
    if (trueIntensity > bound + ACCURACY_CHECK)
    {
      Rprintf("RejectionSampler::simulationStep(): switching rate > bound.\n");
      Rprintf("trueIntensity = %g, bound = %g\n", trueIntensity, bound);
      throw "RejectionSampler::simulationStep(): switching rate > bound.";
    }
    double V = getUniforms(1)(0);
    if (V <= trueIntensity/bound) {
      simulateJump();
      acceptProposedEvent = true;
    }
  }
  else { // i.e. if rejectionFree == true
    simulateJump();
    acceptProposedEvent = true;
  }
    
  updateBound(); // possibly using stored information trueIntensity, proposedEvent, acceptProposedEvent
  return acceptProposedEvent;
}


void ZZAffineRejectionSampler::proposeEvent() {
  
  VectorXd U(getUniforms(dim));
  SizeType index = - 1;
  double deltaT = -1;
  
  for (int i = 0; i < dim; ++i) {
    double simulatedTime = getTimeAffineBound(a(i), b(i), U(i));
    if (simulatedTime > 0 && (index == -1 || simulatedTime < deltaT)) {
      index = i;
      deltaT = simulatedTime;
    }
  }
  if (deltaT < 0)
  {
    throw "ZZAffineRejectionSampler::proposeEvent(): wandered off to infinity.";
  }
  else {
    a += b * deltaT;
    state.x += deltaT * state.v;
    state.t += deltaT;
    proposedEvent = index;
  }
}

Skeleton ZigZag(Sampler& sampler, const int n_iter, const double finalTime) {
  
  sampler.Initialize();
  int iteration = 0;
  Skeleton skel(sampler.getDim(), n_iter);
  skel.Push(sampler.getState());

  while (sampler.getState().t < finalTime || iteration < n_iter) {
    iteration++;
    if(sampler.simulationStep()) // i.e. a switch is accepted
    {
      skel.Push(sampler.getState(), finalTime);
    }
  }
  skel.ShrinkToCurrentSize();
//  Rprintf("ZigZag: Fraction of accepted switches: %g\n", double(skel.getSize()-1)/(iteration));
  return skel;
}


// HELPER FUNCTIONS


double getTimeAffineBound(double a, double b, double u) {
  // simulate T such that P(T>= t) = exp(-at-bt^2/2), using uniform random input u
  // NOTE: Return value -1 indicates +Inf!
  if (b > 0) {
    if (a < 0) 
      return -a/b + getTimeAffineBound(0, b, u);
    else       // a >= 0
      return -a/b + sqrt(a*a/(b * b) - 2 * log(u)/b);
  }
  else if (b == 0) {
    if (a > 0)
      return -log(u)/a;
    else
      return -1; // infinity
  }
  else {
    // b  < 0
    if (a <= 0)
      return -1; // infinity
    else {
      // a > 0
      double t1 = -a/b;
      if (-log(u) <= a * t1 + b * t1 * t1/2)
        return -a/b - sqrt(a*a/(b * b) - 2 * log(u)/b);
      else
        return -1;
    }
  }
}

VectorXd newton(VectorXd& x, const FunctionObject& fn, double precision, const int max_iter) {
  VectorXd grad(fn.gradient(x));
  int i = 0;
  for (i = 0; i < max_iter; ++i) {
    if (grad.norm() < precision)
      break;
    MatrixXd H(fn.hessian(x));
    x -= H.ldlt().solve(grad);
    grad = fn.gradient(x);
  }
  if (i == max_iter) {
    messageStream << "newton: warning: Maximum number of iterations reached without convergence." << std::endl;
  }
  else
    messageStream << "newton: Converged to local minimum in " << i << " iterations." << std::endl;
  return grad;
}

VectorXd resampleVelocity(const int dim, const bool unit_velocity)
{
  // helper function for BPS
  VectorXd v = getStandardNormals(dim);
  if (unit_velocity)
    v.normalize();
  return v;
}
