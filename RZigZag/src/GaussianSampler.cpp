// GaussianSampler.cpp
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

#include "GaussianSampler.h"

bool Gaussian_ZZ::simulationStep() {
  VectorXd U(getUniforms(dim));
  int index = -1;
  double deltaT = -1;
  
  for (int i = 0; i < dim; ++i) {
    double simulatedTime = getTimeAffineBound(a(i), b(i), U(i));
    if (simulatedTime > 0 && (index == -1 || simulatedTime < deltaT)) {
      index = i;
      deltaT = simulatedTime;
    }
  }
  state.x += deltaT * state.v; // O(d)
  state.v[index] = -state.v[index];
  state.t += deltaT;
  z = z + w * deltaT; // O(d), invariant z = V * x
  w = w + 2 * state.v(index) * V.col(index).array(); // preserve invariant w = V * v, O(d)
  a = state.v.array() * z; // invariant a = v .* z, O(d)
  b = state.v.array() * w; // invariant ab = v .* w, O(d)
  
  return true;
}

bool Gaussian_BPS::simulationStep() {
  double t_reflect, t_refresh, deltaT;
  if (refresh_rate <= 0) {
    t_reflect = getTimeAffineBound(a, b, getUniforms(1)(0));
    t_refresh = -1;
    deltaT = t_reflect;
  }
  else {
    VectorXd U(getUniforms(2));
    t_reflect = getTimeAffineBound(a, b, U(0));
    t_refresh = -log(U(1))/refresh_rate;
    deltaT = (t_reflect < t_refresh ? t_reflect : t_refresh);
  }
  state.x += deltaT * state.v; // O(d)
  gradient = gradient + deltaT * w; // O(d)
  state.t += deltaT;
  if (t_refresh < 0 || t_reflect < t_refresh)
  { // O(d)
    VectorXd normalized_gradient = gradient.normalized(); // for projection
    VectorXd delta_v = - 2 * (state.v.dot(normalized_gradient)) * normalized_gradient;
    state.v = state.v + delta_v;
  }
  else
    state.v = resampleVelocity(dim, unit_velocity);
  
  w = V * state.v; // preserves invariant for w, O(d^2)!
  a = state.v.dot(gradient);
  b = state.v.dot(w);

  return true;
}
