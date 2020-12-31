#ifndef __INITIALGUESS_H__
#define __INITIALGUESS_H__

#include <iostream>
#include <math.h>
#include <time.h>

#include <RcppEigen.h>

//[[Rcpp::depends(RcppEigen)]]

//using namespace Eigen;
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

void initialGuess(VectorXd &Theta)
{
  Theta.fill(0.2/Theta.size());
  return;
}

#endif
