// VariousSamplers.h
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

#ifndef __VARIOUSSAMPLERS_H
#define __VARIOUSSAMPLERS_H

#include "ZigZag.h"

class IID_ZZ : public Sampler {
public:
  IID_ZZ(State initialState, double x0 = 0): Sampler(initialState), x0{x0} { };
  IID_ZZ(SizeType dim, double x0 = 0): Sampler(dim), x0{x0} {};
  bool simulationStep();
  void Initialize();
  
  virtual double inversePotentialPlus(double) const = 0;
  virtual double inversePotentialMinus(double) const = 0;
  virtual double univariatePotential(double) const = 0;

protected:
  const double x0;
  double sampleEventTime(const double x, const double v, const double uniform) const;
  VectorXd eventTimes;
};

class Symmetric_IID_ZZ : public IID_ZZ {
public:
  Symmetric_IID_ZZ(State initialState): IID_ZZ(initialState) {};
  Symmetric_IID_ZZ(SizeType dim): IID_ZZ(dim) {};
  double inversePotentialMinus(double y) const  { return -inversePotentialPlus(y);};
};

class StudentT_IID_ZZ : public Symmetric_IID_ZZ {
public:
  StudentT_IID_ZZ(State initialState, double dof): Symmetric_IID_ZZ(initialState), dof{dof} {};
  StudentT_IID_ZZ(SizeType dim, double dof): Symmetric_IID_ZZ(dim), dof{dof}  {};
  
  double univariatePotential(double x) const { return (dof + 1)/2 * log(1 + x*x/dof);}
  double inversePotentialPlus(double y) const { return sqrt(dof * (exp(2 * y/(dof+1)) -1));}
  
private:
  const double dof;
};

class Gaussian_IID_ZZ : public Symmetric_IID_ZZ {
public:
  Gaussian_IID_ZZ(State initialState, double variance): Symmetric_IID_ZZ(initialState), variance{variance} {};
  Gaussian_IID_ZZ(SizeType dim, double variance): Symmetric_IID_ZZ(dim), variance{variance} {};
  
  double univariatePotential(double x) const { return x*x/(2 * variance);}
  double inversePotentialPlus(double y) const { return sqrt(2 * variance * y);}
  
private:
  const double variance;
};


class SphericallySymmetricZZSampler : public ZZAffineRejectionSampler {
  // the setting is U(x) = h(|x|), with h'(z) >= 0 for z >= 0, 
  // and we assume h'(z) \leq C_1  + C_2 z for z >=0
  
public:
  SphericallySymmetricZZSampler(State initialState): ZZAffineRejectionSampler(initialState) { };
  SphericallySymmetricZZSampler(SizeType dim): ZZAffineRejectionSampler(dim) { };
  void Initialize();

protected:
  void updateBound();
  double C1, C2;
  
private:
  double getTrueIntensity();
  virtual double h_prime(double x) const = 0; // 
};

class SphericallySymmetricStudentZZ : public SphericallySymmetricZZSampler {
public:
  // SphericallySymmetricStudentZZ(State initialState, const double dof) : SphericallySymmetricZZSampler(initialState), dof{dof} { C1 = (dof+dim)/(2*sqrt(dof)); C2 = 0;};
  // SphericallySymmetricStudentZZ(SizeType dim, const double dof) : SphericallySymmetricZZSampler(dim), dof{dof} { C1 = (dof+dim)/(2*sqrt(dof)); C2 = 0;};
  SphericallySymmetricStudentZZ(State initialState, const double dof) : SphericallySymmetricZZSampler(initialState), dof{dof} { C1 = 0; C2 = (dof+dim)/dof;};
  SphericallySymmetricStudentZZ(SizeType dim, const double dof) : SphericallySymmetricZZSampler(dim), dof{dof}  { C1 = 0; C2 = (dof+dim)/dof;};
  
private:
  double h_prime(double x) const { return (dof + dim) * x / (dof + x*x);}  

  const double dof;
};

class Affine_BPS : public RejectionSampler {
  // n_clocks = 2, the first of the parallel clocks (i.e. with index 0) represents refreshment
public:
  Affine_BPS(State initialState, const double refresh_rate = 1.0, const bool unit_velocity = false): RejectionSampler(initialState,2), gradient{VectorXd::Zero(initialState.x.size())}, refresh_rate(refresh_rate), unit_velocity{unit_velocity}  { a = VectorXd(2); b = VectorXd(2); a[0] = refresh_rate; b[0] = 0.0;  };
  Affine_BPS(SizeType dim, const double refresh_rate = 1.0, const bool unit_velocity = false): RejectionSampler(dim,2), gradient{VectorXd::Zero(dim)}, refresh_rate(refresh_rate), unit_velocity{unit_velocity}  {  a = VectorXd(2); b = VectorXd(2); a[0] = refresh_rate; b[0] = 0.0;  };

protected:
  double getTrueIntensity();

  VectorXd gradient;
  const double refresh_rate;
  VectorXd a, b;
  
private:
  void proposeEvent();
  double getBound() const { return a[proposedEvent]; }
  virtual void updateGradient() = 0;
  void simulateJump();

  const bool unit_velocity;
};

class IID_BPS : public Affine_BPS {
  // setting is U(x) = sum_i h(x_i), with |h''(z)| \leq M for all z

public:
  IID_BPS(State initialState, const double refresh_rate = 1.0, const bool unit_velocity = false): Affine_BPS(initialState, refresh_rate, unit_velocity) {  };
  IID_BPS(SizeType dim, const double refresh_rate = 1.0, const bool unit_velocity = false): Affine_BPS(dim, refresh_rate, unit_velocity) { } ;
  
protected:
  void updateBound() { a[1] = gradient.dot(state.v); b[1] = state.v.squaredNorm() * M; };
  double M;

private:
  void updateGradient() { for (int i = 0; i < dim; ++i) { gradient[i] = h_prime(state.x(i)); } } // O(d) computation
  void Initialize() { updateBound(); }
  virtual double h_prime(const double z) const = 0; //
};


class StudentT_IID_BPS : public IID_BPS {
public:
  StudentT_IID_BPS(State initialState, double dof, const double refresh_rate = 1.0, const bool unit_velocity = false): IID_BPS(initialState, refresh_rate, unit_velocity), dof{dof} { M = (dof+1)/dof;};
  StudentT_IID_BPS(SizeType dim, double dof, const double refresh_rate = 1.0, const bool unit_velocity = false): IID_BPS(dim, refresh_rate, unit_velocity), dof{dof} { M = (dof+1)/dof;};

private:
  double h_prime(const double z) const { return (dof + 1)*z/(dof + z*z);};
  const double dof;
};


class Gaussian_IID_BPS : public IID_BPS {
public:
  Gaussian_IID_BPS(State initialState, double variance, const double refresh_rate = 1.0, const bool unit_velocity = false): IID_BPS(initialState, refresh_rate, unit_velocity), variance{variance} { M = 1/variance;};
  Gaussian_IID_BPS(SizeType dim, double variance, const double refresh_rate = 1.0, const bool unit_velocity = false): IID_BPS(dim, refresh_rate, unit_velocity), variance{variance} { M = 1/variance;};
  
private:
  void Initialize() { updateBound(); rejectionFree = true;}
  double h_prime(const double z) const { return z/variance;}
  const double variance;
};

class SphericallySymmetricStudentBPS : public Affine_BPS {
// the setting is U(x) = h(|x|), with h'(z) >= 0 for z >= 0, 
// and we assume h'(z) \leq C_1  + C_2 z for z >=0
// here h is -log density of a spherically symmetric student T  with dof degrees of freedom

public:
  SphericallySymmetricStudentBPS(State initialState, const double dof, const double refresh_rate = 1.0, const bool unit_velocity = false) : Affine_BPS(initialState, refresh_rate, unit_velocity), dof{dof}, C1{(dof+dim)/(2*sqrt(dof))}, C2{0} { };
  SphericallySymmetricStudentBPS(SizeType dim, const double dof, const double refresh_rate = 1.0, const bool unit_velocity = false) : Affine_BPS(dim, refresh_rate, unit_velocity), dof{dof}, C1{(dof+dim)/(2*sqrt(dof))}, C2{0}  { };
  // SphericallySymmetricStudentBPS(State initialState, const double dof, const double refresh_rate = 1.0, const bool unit_velocity = false) : Affine_BPS(initialState, refresh_rate, unit_velocity), dof{dof}, C1{0}, C2{(dof+dim)/dof} { };
  // SphericallySymmetricStudentBPS(SizeType dim, const double dof, const double refresh_rate = 1.0, const bool unit_velocity = false) : Affine_BPS(dim, refresh_rate, unit_velocity), dof{dof}, C1{0}, C2{(dof+dim)/dof}  { };

protected:
  void updateBound() { a[1] = state.v.norm() * (C1 + C2 * state.x.norm()); b[1] = state.v.squaredNorm() * C2; };
  
private:
  double h_prime(double z) const { return (dof + dim) * z / (dof + z*z);}  
  void updateGradient() { double x_norm = state.x.norm(); gradient = h_prime(x_norm)/x_norm * state.x; } // O(d) computation
  void Initialize() { updateBound(); }
  
  const double dof, C1, C2;
};

#endif
