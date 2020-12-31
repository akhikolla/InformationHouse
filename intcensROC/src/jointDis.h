#ifndef __JOINTDIS_H__
#define __JOINTDIS_H__

#include <RcppEigen.h>
#include <iostream>
#include "Ispline.h"
#include "Mspline.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

static double jointDis(VectorXd &theta, double t, double m,VectorXd &knot1, VectorXd &knot2 )
{
  int pn = floor(sqrt(theta.size()));
  int qn = floor(sqrt(theta.size()));
  
  double sum = 0;
  
  VectorXd ispline1 = IsplineX(pn, t, knot1);
  VectorXd ispline2 = IsplineX(qn, m, knot2);
  
  for(int i = 0; i < pn; i++)
    for(int j = 0; j < qn; j++)
    {
      sum += theta(i*qn+j)*ispline1(i)*ispline2(j);
    }
    return sum;
}

static double margM(VectorXd &theta, double m, VectorXd &knot2 )
{
  int pn = floor(sqrt(theta.size()));
  int qn = floor(sqrt(theta.size()));
  double alpha_sum = 0;
  double sum = 0;
  VectorXd ispline = IsplineX(qn, m, knot2);
  for(int j = 0; j<pn; j++)
  {
    alpha_sum = 0;
    for(int i = 0; i<qn; i++)
    {
      alpha_sum +=theta(i*qn+j);
    }
    alpha_sum += theta(qn*pn+j);
    sum += alpha_sum*ispline(j);
  }
  return sum;
}

#endif
