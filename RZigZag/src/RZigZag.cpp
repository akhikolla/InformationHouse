// RZigZag.cpp : implements Zig-Zag and other PDMP samplers
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

#include "RZigZag.h"

List SkeletonToList(const Skeleton& skel) {
  // output: R list consisting of Times, Positions and Velocities
  return List::create(Named("Times") = skel.getTimes(), Named("Positions") = skel.getPositions(), Named("Velocities") = skel.getVelocities());
}

Skeleton ListToSkeleton(const List& list) {
  return Skeleton(list["Times"], list["Positions"], list["Velocities"]);
}

//' ZigZagGaussian
//' 
//' Applies the Zig-Zag Sampler to a Gaussian target distribution, as detailed in Bierkens, Fearnhead, Roberts, The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data, 2016.
//' Assume potential of the form \deqn{U(x) = (x - mu)^T V (x - mu)/2,} i.e. a Gaussian with mean vector \code{mu} and covariance matrix \code{inv(V)}
//'
//' @param V the inverse covariance matrix (or precision matrix) of the Gaussian target distribution.
//' @param mu mean of the Gaussian target distribution
//' @param n_iter Number of algorithm iterations; will result in the equivalent amount of skeleton points in Gaussian case because no rejections are needed.
//' @param finalTime If provided and nonnegative, run the sampler until a trajectory of continuous time length finalTime is obtained (ignoring the value of \code{n_iterations})
//' @param x0 starting point (optional, if not specified taken to be the origin)
//' @param v0 starting direction (optional, if not specified taken to be +1 in every component)
//' @return Returns a list with the following objects:
//' @return \code{Times}: Vector of switching times
//' @return \code{Positions}: Matrix whose columns are locations of switches. The number of columns is identical to the length of \code{skeletonTimes}. Be aware that the skeleton points themselves are NOT samples from the target distribution.
//' @return \code{Velocities}: Matrix whose columns are velocities just after switches. The number of columns is identical to the length of \code{skeletonTimes}.
//' @examples
//' V <- matrix(c(3,1,1,3),nrow=2)
//' mu <- c(2,2)
//' result <- ZigZagGaussian(V, mu, 100)
//' plot(result$Positions[1,], result$Positions[2,],type='l',asp=1)
//' @export
// [[Rcpp::export]]
List ZigZagGaussian(const Eigen::MatrixXd V, const Eigen::VectorXd mu, int n_iter = -1, double finalTime = -1, const NumericVector x0 = NumericVector(0), const NumericVector v0 = NumericVector(0)) {
  if (finalTime >= 0)
    n_iter = -1;
  else if (n_iter >= 0)
    finalTime = -1;
  else
    stop("Either finalTime or n_iter must be specified.");
  
  const int dim = V.rows();
  VectorXd x, v;
  if (x0.size() < dim)
    x = VectorXd::Zero(dim);
  else
    x = as<Eigen::Map<VectorXd> >(x0);
  if (v0.size() < dim)
    v = VectorXd::Ones(dim);
  else
    v = as<Eigen::Map<VectorXd> >(v0);
  
  Gaussian_ZZ sampler(V, x, v, mu);
  Skeleton skel(ZigZag(sampler, n_iter, finalTime));
  return SkeletonToList(skel);
}

//' ZigZagLogistic
//'
//' Applies the Zig-Zag Sampler to logistic regression, as detailed in Bierkens, Fearnhead, Roberts, The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data, 2019.
//'
//' @param dataX Design matrix containing observations of the independent variables x. The i-th row represents the i-th observation with components x_{i,1}, ..., x_{i,d}.
//' @param dataY Vector of length n containing {0, 1}-valued observations of the dependent variable y.
//' @param n_iter Number of algorithm iterations; will result in the equivalent amount of skeleton points in Gaussian case because no rejections are needed.
//' @param finalTime If provided and nonnegative, run the sampler until a trajectory of continuous time length finalTime is obtained (ignoring the value of \code{n_iterations})
//' @param x0 starting point (optional, if not specified taken to be the origin)
//' @param v0 starting direction (optional, if not specified taken to be +1 in every component)
//' @param cv optional boolean to indicate the use of subsampling with control variates
//' @return Returns a list with the following objects:
//' @return \code{Times}: Vector of switching times
//' @return \code{Positions}: Matrix whose columns are locations of switches. The number of columns is identical to the length of \code{skeletonTimes}. Be aware that the skeleton points themselves are NOT samples from the target distribution.
//' @return \code{Velocities}: Matrix whose columns are velocities just after switches. The number of columns is identical to the length of \code{skeletonTimes}.
//' @examples
//' require("RZigZag")
//'
//' generate.logistic.data <- function(beta, n.obs) {
//'   dim <- length(beta)
//'   dataX <- cbind(rep(1.0,n.obs), matrix(rnorm((dim -1) * n.obs), ncol = dim -1));
//'   vals <- dataX %*% as.vector(beta)
//'     generateY <- function(p) { rbinom(1, 1, p)}
//'   dataY <- sapply(1/(1 + exp(-vals)), generateY)
//'     return(list(dataX = dataX, dataY = dataY))
//' }
//'
//' beta <- c(1,2)
//' data <- generate.logistic.data(beta, 1000)
//' result <- ZigZagLogistic(data$dataX, data$dataY, 1000)
//' plot(result$Positions[1,], result$Positions[2,],type='l',asp=1)
//' @export
// [[Rcpp::export]]
List ZigZagLogistic(const Eigen::MatrixXd& dataX, const Eigen::VectorXi& dataY, int n_iter = -1, double finalTime = -1, const NumericVector x0 = NumericVector(0), const NumericVector v0 = NumericVector(0), bool cv = false) {

  if (finalTime >= 0)
    n_iter = -1;
  else if (n_iter >= 0)
    finalTime = -1;
  else
    stop("Either finalTime or n_iter must be specified.");
  
  const int dim = dataX.cols();
  VectorXd x, v;
  if (x0.size() < dim)
    x = VectorXd::Zero(dim);
  else
    x = as<Eigen::Map<VectorXd> >(x0);
  if (v0.size() < dim)
    v = VectorXd::Ones(dim);
  else
    v = as<Eigen::Map<VectorXd> >(v0);
  
  if (cv) {
    LogisticCVZZ sampler(dataX, dataY, x, v);
    Skeleton skel(ZigZag(sampler, n_iter, finalTime));
    return SkeletonToList(skel);
  }
  else {
    LogisticZZ sampler(dataX, dataY, x, v);
    Skeleton skel(ZigZag(sampler, n_iter, finalTime));
    return SkeletonToList(skel);
  }
}


//' ZigZagStudentT
//' 
//' Applies the Zig-Zag Sampler to a Student T distribution (IID or spherically symmetric)
//'
//' @param dof scalar indicating degrees of freedom
//' @param dim dimension
//' @param n_iter Number of algorithm iterations; will result in the equivalent amount of skeleton points in Gaussian case because no rejections are needed.
//' @param finalTime If provided and nonnegative, run the sampler until a trajectory of continuous time length finalTime is obtained (ignoring the value of \code{n_iterations})
//' @param x0 starting point (optional, if not specified taken to be the origin)
//' @param v0 starting direction (optional, if not specified taken to be +1 in every component)
//' @param sphericallySymmetric boolean. If false, sample iid Student T distribution, if true (default) sample spherically summetric Student t dsitribution.
//' @return Returns a list with the following objects:
//' @return \code{Times}: Vector of switching times
//' @return \code{Positions}: Matrix whose columns are locations of switches. The number of columns is identical to the length of \code{skeletonTimes}. Be aware that the skeleton points themselves are NOT samples from the target distribution.
//' @return \code{Velocities}: Matrix whose columns are velocities just after switches. The number of columns is identical to the length of \code{skeletonTimes}.
//' @examples
//' dim = 2
//' dof = 4
//' result <- ZigZagStudentT(dof, dim, n_iter=1000, sphericallySymmetric = TRUE)
//' plot(result$Positions[1,], result$Positions[2,],type='l',asp=1)
//' @export
// [[Rcpp::export]]
List ZigZagStudentT(double dof, int dim = 1, int n_iter = -1, double finalTime = -1, const NumericVector x0 = NumericVector(0), const NumericVector v0 = NumericVector(0), bool sphericallySymmetric = true) {
  if (finalTime >= 0)
    n_iter = -1;
  else if (n_iter >= 0)
    finalTime = -1;
  else
    stop("Either finalTime or n_iter must be specified.");
  
  VectorXd x, v;
  if (x0.size() < dim)
    x = VectorXd::Zero(dim);
  else
    x = as<Eigen::Map<VectorXd> >(x0);
  if (v0.size() < dim)
    v = VectorXd::Ones(dim);
  else
    v = as<Eigen::Map<VectorXd> >(v0);
  if (sphericallySymmetric) {
    SphericallySymmetricStudentZZ sampler(State(x,v),dof);
    Skeleton skel = ZigZag(sampler, n_iter, finalTime);
    return SkeletonToList(skel);
  }
  else {
    StudentT_IID_ZZ sampler(State(x,v), dof);
    Skeleton skel = ZigZag(sampler, n_iter, finalTime);
    return SkeletonToList(skel);
  }
}

//' BPSStudentT
//' 
//' Applies the Zig-Zag Sampler to a Student T distribution (IID or spherically symmetric)
//'
//' @param dof scalar indicating degrees of freedom
//' @param dim dimension
//' @param n_iter Number of algorithm iterations; will result in the equivalent amount of skeleton points in Gaussian case because no rejections are needed.
//' @param finalTime If provided and nonnegative, run the sampler until a trajectory of continuous time length finalTime is obtained (ignoring the value of \code{n_iterations})
//' @param x0 starting point (optional, if not specified taken to be the origin)
//' @param v0 starting direction (optional, if not specified taken to be a random vector)
//' @param sphericallySymmetric boolean. If false, sample iid Student T distribution, if true (default) sample spherically summetric Student t dsitribution.
//' @param refresh_rate \code{lambda_refresh}
//' @param unit_velocity TRUE indicates velocities uniform on unit sphere, FALSE (default) indicates standard normal velocities
//' @return Returns a list with the following objects:
//' @return \code{Times}: Vector of switching times
//' @return \code{Positions}: Matrix whose columns are locations of switches. The number of columns is identical to the length of \code{skeletonTimes}. Be aware that the skeleton points themselves are NOT samples from the target distribution.
//' @return \code{Velocities}: Matrix whose columns are velocities just after switches. The number of columns is identical to the length of \code{skeletonTimes}.
//' @examples
//' dim = 2
//' dof = 4
//' result <- BPSStudentT(dof, dim, n_iter=1000,sphericallySymmetric = TRUE)
//' plot(result$Positions[1,], result$Positions[2,],type='l',asp=1)
//' @export
// [[Rcpp::export]]
List BPSStudentT(double dof, int dim = 1, int n_iter = -1, double finalTime = -1, const NumericVector x0 = NumericVector(0), const NumericVector v0 = NumericVector(0), bool sphericallySymmetric = true, const double refresh_rate = 1, const bool unit_velocity = false) {
  if (finalTime >= 0)
    n_iter = -1;
  else if (n_iter >= 0)
    finalTime = -1;
  else
    stop("Either finalTime or n_iter must be specified.");
  
  VectorXd x, v;
  if (x0.size() < dim)
    x = VectorXd::Zero(dim);
  else
    x = as<Eigen::Map<VectorXd> >(x0);
  if (v0.size() < dim) {
    v = as<Eigen::Map<VectorXd> >(rnorm(dim));
    if (unit_velocity)
      v = v/v.norm();
  }
  else
    v = as<Eigen::Map<VectorXd> >(v0);
  if (sphericallySymmetric) {
    SphericallySymmetricStudentBPS sampler(State(x,v), dof, refresh_rate, unit_velocity);
    Skeleton skel = ZigZag(sampler, n_iter, finalTime);
    return SkeletonToList(skel);
  }
  else {
    StudentT_IID_BPS sampler(State(x,v), dof, refresh_rate, unit_velocity);
    Skeleton skel = ZigZag(sampler, n_iter, finalTime);
    return SkeletonToList(skel);
  }
}


//' ZigZagIIDGaussian
//' 
//' Applies the Zig-Zag Sampler to a IID Gaussian distribution
//'
//' @param variance scalar indicating variance
//' @param dim dimension
//' @param n_iter Number of algorithm iterations; will result in the equivalent amount of skeleton points in Gaussian case because no rejections are needed.
//' @param finalTime If provided and nonnegative, run the sampler until a trajectory of continuous time length finalTime is obtained (ignoring the value of \code{n_iterations})
//' @param x0 starting point (optional, if not specified taken to be the origin)
//' @param v0 starting direction (optional, if not specified taken to be +1 in every component)
//' @return Returns a list with the following objects:
//' @return \code{Times}: Vector of switching times
//' @return \code{Positions}: Matrix whose columns are locations of switches. The number of columns is identical to the length of \code{skeletonTimes}. Be aware that the skeleton points themselves are NOT samples from the target distribution.
//' @return \code{Velocities}: Matrix whose columns are velocities just after switches. The number of columns is identical to the length of \code{skeletonTimes}.
//' @examples
//' result <- ZigZagIIDGaussian(1, 2, 1000)
//' plot(result$Positions[2,], result$Positions[1,],type='l',asp=1)
//' @export
// [[Rcpp::export]]
List ZigZagIIDGaussian(double variance, int dim = 1, int n_iter = -1, double finalTime = -1, const NumericVector x0 = NumericVector(0), const NumericVector v0 = NumericVector(0)) {
  if (finalTime >= 0)
    n_iter = -1;
  else if (n_iter >= 0)
    finalTime = -1;
  else
    stop("Either finalTime or n_iter must be specified.");
  
  VectorXd x, v;
  if (x0.size() < dim)
    x = VectorXd::Zero(dim);
  else
    x = as<Eigen::Map<VectorXd> >(x0);
  if (v0.size() < dim)
    v = VectorXd::Ones(dim);
  else
    v = as<Eigen::Map<VectorXd> >(v0);
  
  Gaussian_IID_ZZ sampler(State(x,v), variance);
  Skeleton skel(ZigZag(sampler, n_iter, finalTime));
  return SkeletonToList(skel);
}

//' BPSGaussian
//'
//' Applies the BPS Sampler to a Gaussian target distribution, as detailed in Bouchard-Côté et al, 2017.
//' Assume potential of the form \deqn{U(x) = (x - mu)^T V (x - mu)/2,} i.e. a Gaussian with mean vector \code{mu} and covariance matrix \code{inv(V)}
//'
//' @param V the inverse covariance matrix (or precision matrix) of the Gaussian target distribution.
//' @param mu mean of the Gaussian target distribution
//' @param n_iter Number of algorithm iterations; will result in the equivalent amount of skeleton points in Gaussian case because no rejections are needed.
//' @param x0 starting point (optional, if not specified taken to be the origin)
//' @param v0 starting direction (optional, if not specified taken to be a random vector)
//' @param finalTime If provided and nonnegative, run the BPS sampler until a trajectory of continuous time length finalTime is obtained (ignoring the value of \code{n_iterations})
//' @param refresh_rate \code{lambda_refresh}
//' @param unit_velocity TRUE indicates velocities uniform on unit sphere, FALSE (default) indicates standard normal velocities
//' @return Returns a list with the following objects:
//' @return \code{Times}: Vector of switching times
//' @return \code{Positions}: Matrix whose columns are locations of switches. The number of columns is identical to the length of \code{skeletonTimes}. Be aware that the skeleton points themselves are NOT samples from the target distribution.
//' @return \code{Velocities}: Matrix whose columns are velocities just after switches. The number of columns is identical to the length of \code{skeletonTimes}.
//' @examples
//' V <- matrix(c(3,1,1,3),nrow=2)
//' mu <- c(2,2)
//' x0 <- c(0,0)
//' result <- BPSGaussian(V, mu, n_iter = 100, x0 = x0)
//' plot(result$Positions[1,], result$Positions[2,],type='l',asp=1)
//' @export
// [[Rcpp::export]]
List BPSGaussian(const Eigen::MatrixXd V, const Eigen::VectorXd mu, int n_iter = -1, double finalTime = -1.0, const NumericVector x0 = NumericVector(0), const NumericVector v0 = NumericVector(0), const double refresh_rate = 1, const bool unit_velocity = false) {
  
  if (finalTime >= 0)
    n_iter = -1;
  else if (n_iter >= 0)
    finalTime = -1;
  else
    stop("Either finalTime or n_iterations must be specified.");
  
  const int dim = V.rows();
  VectorXd x, v;
  if (x0.size() < dim)
    x = VectorXd::Zero(dim);
  else
    x = as<Eigen::Map<VectorXd> >(x0);
  if (v0.size() < dim) {
    v = as<Eigen::Map<VectorXd> >(rnorm(dim));
    if (unit_velocity)
      v = v/v.norm();
  }
  else
    v = as<Eigen::Map<VectorXd> >(v0);
  
  Gaussian_BPS sampler(V, x, v, mu, refresh_rate, unit_velocity);
  Skeleton skel(ZigZag(sampler, n_iter, finalTime));
  return SkeletonToList(skel);
}

//' BPSIIDGaussian
//' 
//' Applies the Bouncy Particle Sampler to a IID Gaussian distribution
//'
//' @param variance scalar indicating variance
//' @param dim dimension
//' @param n_iter Number of algorithm iterations; will result in the equivalent amount of skeleton points in Gaussian case because no rejections are needed.
//' @param finalTime If provided and nonnegative, run the sampler until a trajectory of continuous time length finalTime is obtained (ignoring the value of \code{n_iterations})
//' @param x0 starting point (optional, if not specified taken to be the origin)
//' @param v0 starting direction (optional, if not specified taken to be a random vector)
//' @param refresh_rate \code{lambda_refresh}
//' @param unit_velocity TRUE indicates velocities uniform on unit sphere, FALSE (default) indicates standard normal velocities
//' @return Returns a list with the following objects:
//' @return \code{Times}: Vector of switching times
//' @return \code{Positions}: Matrix whose columns are locations of switches. The number of columns is identical to the length of \code{skeletonTimes}. Be aware that the skeleton points themselves are NOT samples from the target distribution.
//' @return \code{Velocities}: Matrix whose columns are velocities just after switches. The number of columns is identical to the length of \code{skeletonTimes}.
//' @examples
//' result <- BPSIIDGaussian(1, 2, 1000)
//' plot(result$Positions[2,], result$Positions[1,],type='l',asp=1)
//' @export
// [[Rcpp::export]]
List BPSIIDGaussian(double variance, int dim = 1, int n_iter = -1, double finalTime = -1, const NumericVector x0 = NumericVector(0), const NumericVector v0 = NumericVector(0), const double refresh_rate = 1, const bool unit_velocity = false) {
  if (finalTime >= 0)
    n_iter = -1;
  else if (n_iter >= 0)
    finalTime = -1;
  else
    stop("Either finalTime or n_iter must be specified.");
  
  VectorXd x, v;
  if (x0.size() < dim)
    x = VectorXd::Zero(dim);
  else
    x = as<Eigen::Map<VectorXd> >(x0);
  if (v0.size() < dim) {
    v = as<Eigen::Map<VectorXd> >(rnorm(dim));
    if (unit_velocity)
      v = v/v.norm();
  }
  else
    v = as<Eigen::Map<VectorXd> >(v0);
  
  Gaussian_IID_BPS sampler(State(x,v), variance, refresh_rate, unit_velocity);
  Skeleton skel(ZigZag(sampler, n_iter, finalTime));
  return SkeletonToList(skel);
}

//' EstimateESS
//' 
//' Estimates the effective sample size (ESS) of a piecewise deterministic skeleton
//' 
//' @param skeletonList a piecewise deterministic skeleton (consisting of Times, Points and Velocities) returned by a sampler
//' @param n_batches optional argument indicating the number of batches to use in the batch means estimation method
//' @param coordinate if specified, only estimate the ESS of the specified coordinate, otherwise estimate the ESS of all coordinates
//' @param zeroMeans if TRUE do not estimate means but assume a centered distribution
//' @return Returns a list containing the estimated asymptotic variance, ESS and estimated covariance matrix
//' @export
// [[Rcpp::export]]
List EstimateESS(const List& skeletonList, int n_batches = 100, int coordinate = -1, bool zeroMeans = false) {
  
  Skeleton skel = ListToSkeleton(skeletonList);
  PostProcess pp(skel);
  if (coordinate > 0)
    coordinate--; // convert R to C++ number
  pp.estimateESS(n_batches, coordinate, zeroMeans);
  return List::create(Named("AsVar") = pp.getAsVar(), Named("ESS") = pp.getESS(), Named("Cov") = pp.getCovarianceMatrix());
}

//' EstimateMoment
//' 
//' Estimates the p-th moment of a piecewise deterministic skeleton
//' 
//' @param skeletonList a piecewise deterministic skeleton (consisting of Times, Points and Velocities) returned by a sampler
//' @param p moment to estimate
//' @param coordinate if specified, only estimate the ESS of the specified coordinate, otherwise estimate the ESS of all coordinates
//' @return Returns a list containing the estimated moment
//' @export
// [[Rcpp::export]]
List EstimateMoment(const List& skeletonList, const int p, int coordinate = -1) {
  
  Skeleton skel = ListToSkeleton(skeletonList);
  PostProcess pp(skel);
  if (coordinate > 0)
    coordinate--; // convert R to C++ number
  pp.estimateMoment(p, coordinate);
  return List::create(Named("Moment") = pp.getMoment());
}


//' EstimateCovarianceMatrix
//' 
//' Estimates the covariance matrix of a piecewise deterministic skeleton
//' 
//' @param skeletonList a piecewise deterministic skeleton (consisting of Times, Points and Velocities) returned by a sampler
//' @param coordinate if specified, only estimate the variance of the specified coordinate, otherwise estimate the covariance matrix of all coordinates
//' @param zeroMeans if TRUE do not estimate means but assume a centered distribution
//' @return Returns a list containing the estimated moment
//' @export
// [[Rcpp::export]]
List EstimateCovarianceMatrix(const List& skeletonList, int coordinate = -1, bool zeroMeans = false) {
  
  Skeleton skel = ListToSkeleton(skeletonList);
  PostProcess pp(skel);
  if (coordinate > 0)
    coordinate--; // convert R to C++ number
  pp.estimateCovariance(coordinate, zeroMeans);
  return List::create(Named("Cov") = pp.getCovarianceMatrix());
}

//' DiscreteSamples
//' 
//' Extract discrete samples from a skeleton
//' 
//' @param skeletonList a piecewise deterministic skeleton (consisting of Times, Points and Velocities) returned by a sampler
//' @param n_samples number of samples to obtain
//' @param coordinate if specified, only obtain samples of the specified coordinate, otherwise obtain samples of all coordinates
//' @return Returns a list containing the extracted samples and the times (on the continuous time scale) at which the samples are extracted
//' @export
// [[Rcpp::export]]
List DiscreteSamples(const List& skeletonList, const int n_samples, int coordinate = -1) {
  
  Skeleton skel = ListToSkeleton(skeletonList);
  PostProcess pp(skel);
  if (coordinate > 0)
    coordinate--; // convert R to C++ number
  pp.sample(n_samples, coordinate);
  return List::create(Named("Samples") = pp.getSamples(),Named("Times") = pp.getSampleTimes());
}
