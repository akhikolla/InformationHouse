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
double falsePos(NumericVector thetaTmp, double m, double max_m, double t,NumericVector KnotI, NumericVector KnotM)
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
  double jointProMax_m = 0.0;
  double margPro_m = 0.0;
  jointPro = jointDis(theta, t, m, knot1, knot2);
  jointProMax_m = jointDis(theta, t, max_m, knot1, knot2);
  margPro_m = margM(theta,m, knot2);
  return (margPro_m-jointPro)/(1-jointProMax_m);
}
