/*function to calculate true positive rate*/
#include <iostream>
#include "jointDis.h"
#include <Rcpp.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
double truePos(NumericVector thetaTmp, double mtmp, double max_m_tmp, double ttmp,NumericVector KnotI, NumericVector KnotM)
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
  
  double t = ttmp;
  double m = mtmp;
  double max_m = max_m_tmp; 
  
  
  double jointPro = 0.0;
  double jointProMax_m = 0.0;
  jointPro = jointDis(theta, t, m, knot1, knot2);
  jointProMax_m = jointDis(theta, t, max_m, knot1, knot2);
  return jointPro/jointProMax_m;
}

