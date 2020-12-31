/*
 Jiaxing Lin
 Aug 25 2016
 */

/*
 subroutine to perform the updating of the constraint factors gamma.
 */

#include "global.h"

using namespace Eigen;
using namespace std;

double conGamma(VectorXd &fdirection, VectorXd &theta2)
{
  double threadHolda = 1e-8;
  double threadHoldc = 1e-5;
  double threadHoldb = 0;
  int size = theta2.size(); // number of variables
  double threadHold_coverge = 1e-8;
  // thread hold to check the stability of gamma 
	double gamma_threadhold = 1e-8;	
	double gamma = 1e100;
  
  int  col = floor(sqrt(size));
  
	int  row = size/col;
  row = row + 1;
  
  for(int i = 0; i < col; i++)
  {
    if(fdirection(i) < 0 && theta2(i) > threadHoldb)
      if(gamma > -theta2(i)/fdirection(i))
      {
        gamma = -theta2(i)/fdirection(i);
        *min_lamba = i;
      }
  }
  
  for(int i = 1; i < col; i++)
  {
    if(fdirection(i*col) < 0 && theta2(i*col) > threadHoldb)
      if(gamma > -theta2(i*col)/fdirection(i*col))
      {
        gamma = -theta2(i*col)/fdirection(i*col);
        *min_lamba = i*col;
      }
  }
  
  for(int i = 1; i < col; i++)
  {
    for(int j = 1; j < col; j++)
    {
      if(fdirection(i*col+j) < 0 && theta2(i*col+j) > threadHoldc)
        if(gamma > (threadHoldc-theta2(i*col+j))/fdirection(i*col+j))
        {
          gamma = (threadHoldc-theta2(i*col+j))/fdirection(i*col+j);
          *min_lamba =i*col+j;
          *cse = 1;
        }
    }
  }
  
  for(int i = size-col+1; i < size; i++)
  {
    if(fdirection(i) < 0 && theta2(i) > threadHoldb)
      gamma = min(gamma, -theta2(i)/fdirection(i));
  }
  
  if(fdirection.sum() > threadHolda && (1-theta2.sum()) > threadHolda)
  {  
	  gamma = min(gamma, (1- theta2.sum())/fdirection.sum());
    *min_lamba = -1;
	  *cse = 2;	
  }
  
	if (gamma < 0 && gamma_threadhold < 1e-6)
    gamma = 1;

	
  return min(gamma, 0.5);	
}










