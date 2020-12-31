// ZigZag.h
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

#ifndef __ZIGZAG_H
#define __ZIGZAG_H

#include "RInterface.h"

#ifndef __INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

#define EIGEN_NO_MALLOC

// the following declarations should be defined in an interface file, such as RInterface.h
// this should be the only adapation if we are not compiling using Rcpp
#ifndef __RINTERFACE_H
VectorXd getUniforms(const long n);
VectorXd getStandardNormals(const long n);
#include <iostream>
#include <cmath>
static std::ostream& messageStream = std::cout;
#endif

#define DEFAULTSIZE 1e4

typedef Eigen::Index SizeType;

struct State {
  State(const SizeType dim) : t{0.0}, x{VectorXd::Zero(dim)}, v{VectorXd::Ones(dim)} {};
  State(double t, const VectorXd& x, const VectorXd& v): t{t}, x{x}, v{v} {};
  State(const VectorXd& x, const VectorXd& v): State(0,x,v) {};
  double t;
  VectorXd x;
  VectorXd v; // represents the velocity exactly at time t
};

class Sampler { // abstract class
public:
  // philosophy: Sampler contains the necessary information about the target distribution and keeps track of the current state of the piecewise deterministic sampler
  Sampler(State initialState): dim{initialState.x.size()}, state{initialState} {};
  Sampler(SizeType dim): dim{dim}, state{State(dim)} {};
  virtual bool simulationStep() = 0;
  const State& getState() const { return state;};
  SizeType getDim() const { return dim;};
  virtual void Initialize() {}; // automatically called in the beginning of the main zigzag routine
  
protected:
  const SizeType dim; 
  State state;
};

class Skeleton {
public:
  Skeleton(const int dim, int initialSize = DEFAULTSIZE);
  Skeleton(const VectorXd& Times, const MatrixXd& Positions, const MatrixXd& Velocities): Times{Times}, Positions{Positions}, Velocities{Velocities}, size{Times.size()}, capacity{size}, dimension{Positions.rows()} {};
  // const State getState(int i) { return State(Times(i), Positions.col(i), Velocities.col(i));};
  // Event getState(const double t) const; // get state at continuous time t, TO BE IMPLEMENTED
  const int getSize() const { return size; };
  const MatrixXd& getPositions() const { return Positions; };
  const MatrixXd& getVelocities() const { return Velocities; };
  const VectorXd& getTimes() const { return Times; };
  
  friend class PostProcess;

private:
  void Push(const State& state, const double finalTime = -1);
  void ShrinkToCurrentSize(); // shrinks to actual size;
  void Resize(const int factor = 2);

  VectorXd Times;
  MatrixXd Positions;
  MatrixXd Velocities;
  SizeType size;
  SizeType capacity;
  SizeType dimension;
  
  friend Skeleton ZigZag(Sampler& sampler, const int n_iter, const double finalTime);
};

// the main zigzag routine
Skeleton ZigZag(Sampler& sampler, const int n_iter, const double finalTime);

class RejectionSampler : public Sampler {
public:
  RejectionSampler(const SizeType dim, const SizeType n_clocks): Sampler(dim), n_clocks{n_clocks}, rejectionFree{false} {};
  RejectionSampler(State initialState, const SizeType n_clocks): Sampler(initialState), n_clocks{n_clocks}, rejectionFree{false} {};
  bool simulationStep();
  
protected:
  virtual void updateBound() = 0; // 
  // the values trueIntensity, proposedEvent and acceptProposedEvent concern the most recently proposedEvent
  double trueIntensity; // this value is stored in the class because it is used for acceptance and possibly updating bound
  SizeType proposedEvent; 
  bool acceptProposedEvent;
  const SizeType n_clocks;
  bool rejectionFree;
  
private:
  virtual void proposeEvent() = 0; // moves the state of the sampler to a proposed event and returns the proposed switch
  virtual double getBound() const = 0;
  virtual double getTrueIntensity() = 0;
  virtual void simulateJump() = 0;
};

class ZZAffineRejectionSampler : public RejectionSampler {
public:
  ZZAffineRejectionSampler(const SizeType dim): RejectionSampler(dim,dim) {};
  ZZAffineRejectionSampler(State initialState): RejectionSampler(initialState,initialState.x.size()) {};

protected:
  ArrayXd a, b;

private:
  void proposeEvent(); 
  void simulateJump() { state.v(proposedEvent) = - state.v(proposedEvent);};
  double getBound() const { return a(proposedEvent);};
};

// various helper functions for zigzag sampling

double getTimeAffineBound(double a, double b, double u);

inline double pospart(const double a) {
  if (a > 0)
    return a;
  else
    return 0;
}

class FunctionObject {
public:
  virtual VectorXd gradient(const VectorXd&) const = 0;
  virtual MatrixXd hessian(const VectorXd&) const = 0;
};

VectorXd newton(VectorXd& x, const FunctionObject& fn, double precision = 1e-10, const int max_iter = 1e2);

VectorXd resampleVelocity(const int dim, const bool unit_velocity);

#endif
