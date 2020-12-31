/*
 Jiaxing Lin
 Aug. 25 2016
 */

/*
 function to compute the gradient for the loglikelihood function
 */
#ifndef __sumterms_h__
#define __sumterms_h__

//#include <Eigen/Core>
#include <iostream>

//using namespace Eigen;
using namespace std;

inline double sumterms1(VectorXd &theta, MatrixXd &ispline_U, MatrixXd &mspline_m, int k)
{
  int pn = ispline_U.col(0).size();
  int qn = mspline_m.col(0).size();
  
  double sum = 0;
  
  for(int i = 0; i < pn; i++)
    for(int j = 0; j < qn; j++)
    {
      sum += theta(i*qn+j)*ispline_U(i,k)*mspline_m(j,k);
    }
    return sum;
}



inline double sumterms2(VectorXd &theta, MatrixXd &ispline_U, MatrixXd &mspline_m, int k)
{
  int pn = ispline_U.col(0).size();	
  int qn = mspline_m.col(0).size();
  double alpha_sum = 0;	
  double sum = 0;
  for(int j = 0; j<pn; j++)
  {
    alpha_sum = 0;
    for(int i = 0; i<qn; i++)
    {	
      alpha_sum +=theta(i*qn+j); 					
    }			
    alpha_sum += theta(qn*pn+j);
    sum += alpha_sum*mspline_m(j,k);
  }
  
  return sum;
}

#endif
