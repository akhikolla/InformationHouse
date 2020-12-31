/* function to calculating the false positive rate */


#include <RcppEigen.h>
#include <iostream>
#include "jointDis.h"
#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;
using namespace Rcpp;


//[[Rcpp::export]]
double surva(NumericVector thetaTmp, double m, double m2, double t,NumericVector KnotI, NumericVector KnotM)
{
  VectorXd theta(thetaTmp.size());
  for(int i = 0; i < thetaTmp.size(); i++)
    theta(i) = thetaTmp(i);
  VectorXd knot1(KnotI.size());
  for(int i = 0; i < KnotI.size(); i++)
    knot1(i) = KnotI(i);
  VectorXd knot2(KnotM.size());
  for(int i = 0; i < KnotM.size(); i++)
    knot2(i) = KnotM(i);

  double jointPro = 0.0;
  double jointPro1 = 0.0;
  double margPro_m = 0.0;
  double margPro_m1 = 0.0;
  jointPro = jointDis(theta, t, m, knot1, knot2);
  jointPro1 = jointDis(theta, t, m2, knot1, knot2);
  margPro_m = margM(theta,m, knot2);
  margPro_m1 = margM(theta,m2, knot2);
  return (jointPro - jointPro1)/(margPro_m - margPro_m1);
}
