/*
 Jiaxing Lin
 Aug. 25 2016
 */

/*
 function to compute the loglikelihood function, gradient, hessian 
 matrix.
 */

#ifndef __likelihood_h__
#define __likelihood_h__

//#include <Eigen/Core>
#include <iostream>
#include <algorithm>
#include <math.h>

#include "sumterms.h"

using namespace Eigen;
using namespace std;

double loglikelihood(VectorXd &theta, MatrixXd &ispline_U, 
       MatrixXd &ispline_V, MatrixXd &mspline_m, VectorXd &delta, int pn, int qn)
{
  // Need to check size of all these matrix.
  double Loglikelihood = 1;
  
  for(int k = 0; k < delta.size(); k++)
  {
    if(delta(k) == 1)
      Loglikelihood += log(sumterms1(theta,ispline_U,mspline_m,k));	
    if(delta(k) == 2)
      Loglikelihood += log(sumterms1(theta,ispline_V, mspline_m, k)-
        sumterms1(theta,ispline_U,mspline_m, k));
    if(delta(k) == 3)		
      Loglikelihood += log(sumterms2(theta,ispline_U,mspline_m,k)
                             -sumterms1(theta,ispline_V,mspline_m,k));
  }
  
  return Loglikelihood;	
}

#endif
