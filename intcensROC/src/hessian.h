/*
 Jiaxing Lin
 Aug. 25 2016
 */

/*
 function to compute the gradient for the loglikelihood function
 */

#ifndef __HESSIAN_H__
#define __HESSIAN_H__

//#include <Eigen/Core>
#include <iostream>
#include <algorithm>
// for tesint of computing logrithm only
#include <math.h>

#include "sumterms.h"

using namespace Eigen;
using namespace std;


MatrixXd hessian(VectorXd &theta, MatrixXd &ispline_U, MatrixXd &ispline_V, MatrixXd &mspline_m, VectorXd &delta, int pn, int qn)
{
  // Need to check size of all these matrix.
  MatrixXd Hessian_mat(theta.size(), theta.size());	
  Hessian_mat.fill(0);
  MatrixXd Check(theta.size(), theta.size());
  Check.fill(0);
  // precompute the sum terms for the donominator to save computional cost 
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
  
  // precompute dorminator terms to save computional cost 
  VectorXd term1(delta.size());
  VectorXd term2(delta.size());
  VectorXd term3(delta.size());
  
  for(int k = 0; k < delta.size(); k++)
  {
    if(delta(k) == 1)
      term1(k) = 1/(sumterms1_U_k(k)*sumterms1_U_k(k));
    else
      term1(k) = 0;
	  
    if(delta(k) == 2)
      term2(k) = 1/((sumterms1_V_k(k)-sumterms1_U_k(k))*(sumterms1_V_k(k)-sumterms1_U_k(k)));
    else
   	  term2(k) = 0;
	
	if(delta(k) == 3)
	  term3(k) = 1/((sumterms2_m(k)-sumterms1_V_k(k))*(sumterms2_m(k)-sumterms1_V_k(k))); 
	else
      term3(k) = 0;
  }

	
  // compute the hessian elements for the h(alphas, alphas)
  // there could be double computation her.
  for(int i0=0; i0< pn; i0++){
    for(int j0=0; j0< qn; j0++){
      for(int i = 0; i < pn; i++){
        for(int j =0; j < qn; j++){	
		  if(Check(i0*qn+j0,i*qn+j) == 0)
		  {
		    for(int k = 0; k < delta.size(); k++)		
            {
              if(delta(k) == 1)
                Hessian_mat(i0*qn+j0, i*qn+j) +=-ispline_U(i0,k)*mspline_m(j0,k)
                                              *ispline_U(i,k)*mspline_m(j,k) *term1(k);
              if(delta(k) == 2)
                Hessian_mat(i0*qn+j0,i*qn+j) += -( ispline_V(i0,k)*mspline_m(j0,k)
                                                  -ispline_U(i0,k)*mspline_m(j0,k))
                								 *(ispline_V(i,k)*mspline_m(j,k)
                    							 - ispline_U(i,k)*mspline_m(j,k))*term2(k);
              if(delta(k) == 3) 
                Hessian_mat(i0*qn+j0,i*qn+j) += -(mspline_m(j0,k) 
                                                -ispline_V(i0,k)*mspline_m(j0,k))
                                                *(mspline_m(j,k) 
                                                - ispline_V(i,k)*mspline_m(j,k))*term3(k);
            }
            Check(i0*qn+j0,i*qn+j) = 1;
		    if(Check(i*qn+j, i0*qn+j0) == 0)
			{
              Hessian_mat(i*qn+j,i0*qn+j0)=Hessian_mat(i0*qn+j0, i*qn+j);
			  Check(i*qn+j,i0*qn+j0) = 1;
            }
          }
        }
      }
	}
  }
  
  // compute the hessian matrix element h(alphas, gammas) and h(gammas, alphas)
  for(int i = 0; i < pn; i++){
    for(int j = 0; j < qn; j++){
      for(int j0 = 0; j0 <qn; j0++){
        for(int k = 0; k < delta.size(); k++){
          if(delta(k) == 3)
            Hessian_mat(i*qn+j,pn*qn+j0) += -mspline_m(j0,k)*(mspline_m(j,k) -		
              ispline_V(i,k)*mspline_m(j,k))*term3(k);
        }
        Hessian_mat(pn*qn+j0, i*qn+j) = Hessian_mat(i*qn+j,pn*qn+j0);
      }
    }
  }
  
  
  // compute the hessian matrix element 
  
  for(int j = 0; j < qn; j++){
    for(int j0 = j; j0 < qn; j0++){
      for(int k = 0; k < delta.size(); k++){
        if(delta(k) == 3)
          Hessian_mat(pn*qn+j,pn*qn+j0) += -mspline_m(j,k)*mspline_m(j0,k)*term3(k);
      }
      Hessian_mat(pn*qn+j0,pn*qn+j) = Hessian_mat(pn*qn+j,pn*qn+j0);
    }
  }
  
  return Hessian_mat;
}


#endif
