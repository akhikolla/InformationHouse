// GaussianSampler.h
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

#ifndef __GAUSSIANSAMPLER_H
#define __GAUSSIANSAMPLER_H
#include "ZigZag.h"

class Gaussian_ZZ : public Sampler {
public:
  Gaussian_ZZ(const MatrixXd& V, const VectorXd& mu = VectorXd(0)): Sampler(State(V.rows())), V{V}, mu{mu} { };
  Gaussian_ZZ(const MatrixXd& V, VectorXd x, VectorXd v, const VectorXd& mu = VectorXd(0)): Sampler(State(x,v)), V{V}, mu{mu} { };
  bool simulationStep();
  void Initialize() {w = (V * state.v).array(); z =(V * (state.x - mu)).array(); a = state.v.array() * z; b = state.v.array() * w;};
  
private:
  const MatrixXd& V;
  const VectorXd& mu;
  ArrayXd w, z; // invariants w = V * v, z = V * x
  ArrayXd a, b; // invariants a = v.* z, b = v .* w
};

class Gaussian_BPS : public Sampler {
public:
  Gaussian_BPS(const MatrixXd& V, const VectorXd& mu = VectorXd(0), const double refresh_rate = 1.0, const bool unit_velocity = false): Sampler(State(V.rows())), V{V}, mu{mu}, refresh_rate(refresh_rate), unit_velocity{unit_velocity} {  };
  Gaussian_BPS(const MatrixXd& V, VectorXd x, VectorXd v, const VectorXd& mu = VectorXd(0), const double refresh_rate = 1.0, const bool unit_velocity = false): Sampler(State(x,v)), V{V}, mu{mu}, refresh_rate(refresh_rate), unit_velocity{unit_velocity} { };
  bool simulationStep();
  void Initialize() { gradient = V * (state.x - mu); w = V * state.v; a = state.v.dot(gradient); b = state.v.dot(w); };
  
private:
  const MatrixXd& V;
  const VectorXd& mu;
  const double refresh_rate;
  const bool unit_velocity;
  VectorXd gradient, w;
  double a, b;
};



#endif
