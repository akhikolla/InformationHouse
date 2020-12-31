/*
 Jiaxing Lin
 Aug. 25 2016
 */

/*
 updating the solutions, that is move further the variables towards
 the maximum of the logrithm likelihood.
 */

#include <iostream>
#include <algorithm>
#include <math.h>

#include "likelihood.h"
#include "global.h"

using namespace Eigen;
using namespace std;

void updateTheta(VectorXd &fdireaction, double gamma, VectorXd &theta,
                 MatrixXd &ispline_U, MatrixXd &ispline_V, MatrixXd &mspline_m, VectorXd &delta,
                 int pn, int qn)
{
  double threadHold = 1e-8;
  double threadHoldc = 1e-5;
  int k_limit = 20; // what if thie k needed is larger that this?
  int k = 0;
  int size = theta.size();
  int col = floor(sqrt(size));
  VectorXd thetaUpdate = theta;
	int row = size/col; 
  row = row + 2;

  double half_pow_k = 1.0;
  
  VectorXd tempTheta;
  // checking for the smalles k that match the updating critiea
  for(int i = 0; i <= k_limit; i++)
  {
    tempTheta = theta + half_pow_k*gamma*fdireaction;
    thetaUpdate = tempTheta * threadHold;  
    if(loglikelihood(tempTheta, ispline_U, ispline_V, mspline_m, delta, pn, qn) 
            >=loglikelihood(theta, ispline_U, ispline_V, mspline_m, delta, pn, qn))
      break;		
    else
    {
      half_pow_k = half_pow_k*0.5;	
      k++;
    }
  }
  
  //update theta
  theta = theta + min(half_pow_k*gamma,0.5) *fdireaction;
  
  for(int i = 0; i < col; i++)
  {
    if(theta(i) < threadHold)
      theta(i) = 0;
    if(theta(i*col) < threadHold)
  	  theta(i*col) = 0;
	  if(theta(size-col+i) < threadHold)
	    theta(size-col+i) = 0;
  }

  for(int i = 1; i < col; i++)
    for(int j = 1; j < col; j++)
      if(abs(theta(i*col+j)-threadHoldc) < threadHold/1000 )
       	theta(i*col+j) = threadHoldc;

  thetaUpdate = tempTheta * threadHold;  
  return;
}








