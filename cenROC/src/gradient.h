/*
 Jiaxing Lin
 Aug. 25 2016
 */

/*
 function to compute the gradient for the loglikelihood function
 */
#ifndef __gradient_h__
#define __gradient_h__



//#include <Eigen/Core>
#include <iostream>
#include <algorithm>
// for tesint of computing logrithm only
#include <math.h>

#include "sumterms.h"

using namespace Eigen;
using namespace std;

VectorXd gradient(VectorXd &theta, MatrixXd &ispline_U, MatrixXd &ispline_V, MatrixXd &mspline_m, VectorXd &delta, int pn, int qn)
{
  double tol = 1e-10; // set up a numerical tol
  // Need to check size of all these matrix.
  VectorXd Gradient(theta.size());	
  Gradient.fill(0); 
  VectorXd sumterms1_U_k(delta.size());
  VectorXd sumterms1_V_k(delta.size());
  VectorXd sumterms2_m(delta.size());
  
  for(int k = 0; k < delta.size(); k++)
  {
    if(delta(k) == 3)
      sumterms1_U_k(k) = 0;
    else
      sumterms1_U_k(k) = sumterms1(theta,ispline_U,mspline_m,k);
    if(delta(k) == 1)
      sumterms1_V_k(k) = 0;
    else
      sumterms1_V_k(k) = sumterms1(theta,ispline_V,mspline_m,k);
    
    sumterms2_m(k) = sumterms2(theta,ispline_U,mspline_m,k);
  }
  
  // There is some problem with the negative values for this logrithm computation
  // The reason is that: there is no gurantee in the value of the values that is 
  // been taken logrithm.
  // will leave this to later on.
  
  // compute the gradient elements for the alphas
  
  for(int i =0; i < pn; i++){
    for(int j =0; j < qn; j++){	
      for(int k = 0; k < delta.size(); k++)		
      {
        if(delta(k) == 1)
          Gradient(i*qn+j) += ispline_U(i,k)*mspline_m(j,k)
          /sumterms1_U_k(k);
        if(delta(k) == 2)
          Gradient(i*qn+j) += (ispline_V(i,k)*mspline_m(j,k)-ispline_U(i,k)*mspline_m(j,k))/
        ( 
            sumterms1_V_k(k)-
              sumterms1_U_k(k)
        );
        if(delta(k) == 3) 
          Gradient(i*qn+j) += (mspline_m(j,k) - ispline_V(i,k)*mspline_m(j,k))/
        (
            sumterms2_m(k)
          -sumterms1_V_k(k)
        );
      }
    }
  }
  
  // compute the gradient elements for gammas
  for(int j =0; j < qn; j++)
  {
    for(int k = 0; k < delta.size(); k++)
    {
      if(delta(k) == 3)
        Gradient(pn*qn+j) += mspline_m(j,k)/
      (
          sumterms2_m(k)-
            sumterms1_V_k(k)
      );
    }
  }
  
 for(int i = 0;  i < Gradient.size(); i++)
    if(Gradient(i) < tol)
      Gradient(i) = 0;
    
  return Gradient;
}


#endif
