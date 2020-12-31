// VariousSamplers.cpp
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

#include "VariousSamplers.h"

void IID_ZZ::Initialize() {
  VectorXd uniforms = getUniforms(dim);
  eventTimes = VectorXd(dim);
  for (int i = 0; i < dim; ++i) {
    eventTimes(i) = sampleEventTime(state.x(i),state.v(i), uniforms(i));
  }
}

double IID_ZZ::sampleEventTime(const double x, const double v, const double uniform) const {
  
  double U_val = (v *(x-x0) > 0 ? univariatePotential(x) : univariatePotential(x0));
  if (v > 0) {
    return -x/v + inversePotentialPlus(U_val - log(uniform))/v;
  }
  else {
    return -x/v + inversePotentialMinus(U_val - log(uniform))/v;
  }
}

bool IID_ZZ::simulationStep() {
  SizeType index;
  double deltaT = eventTimes.minCoeff(&index);
  eventTimes = eventTimes.array() - deltaT;
  state.x += state.v * deltaT;
  state.t += deltaT;
  state.v(index) = -state.v(index);
  eventTimes(index) = sampleEventTime(state.x(index), state.v(index), getUniforms(1)(0));
  return true;
}

double SphericallySymmetricZZSampler::getTrueIntensity() {
  double x_norm = state.x.norm();
  return h_prime(x_norm)/x_norm * pospart(state.x(proposedEvent) * state.v(proposedEvent));
}

void SphericallySymmetricZZSampler::updateBound() {
  a = C2 * state.v.array() * state.x.array() + C1; // O(d) computation

}

void SphericallySymmetricZZSampler::Initialize() {
  updateBound();
  b = C2 * ArrayXd::Ones(dim);
}


void Affine_BPS::proposeEvent() {
  
  VectorXd U(getUniforms(n_clocks));
  SizeType index = - 1;
  double deltaT = -1;
  
  for (int i = 0; i < n_clocks; ++i) {
    double simulatedTime = getTimeAffineBound(a(i), b(i), U(i));
    if (simulatedTime > 0 && (index == -1 || simulatedTime < deltaT)) {
      index = i;
      deltaT = simulatedTime;
    }
  }
  if (deltaT < 0)
  { 
    // Rprintf("Wandered off to infinity\n");
    throw "Affine_BPS::proposeEvent(): wandered off to infinity.";
  }
  else {
    a[1] += b[1] * deltaT;
    state.x += deltaT * state.v;
    state.t += deltaT;
    proposedEvent = index;
  }
}

void Affine_BPS::simulateJump() {
  
  if (proposedEvent == 0) // refreshment 
  {
    state.v = resampleVelocity(dim, unit_velocity);
  }
  else
  { // O(d)
    VectorXd normalized_gradient = gradient.normalized(); // for projection
    VectorXd delta_v = - 2 * (state.v.dot(normalized_gradient)) * normalized_gradient;
    state.v = state.v + delta_v;
  }
}

double Affine_BPS::getTrueIntensity() {
  updateGradient();
  if (proposedEvent == 0) {
    return refresh_rate;
  }
  else {
    return gradient.dot(state.v);
  }
}
