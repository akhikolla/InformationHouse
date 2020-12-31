/*
 Jiaxing Lin
 Aug. 25 2016
 */

/*
 subroutine to updating the matrix X and Alpha array.
 Alpha and X are mapped to each other, only focus on
 one of them is enough. This subroutine updates X
 */

#include <iostream>
#include <unordered_set>

using namespace Eigen;
using namespace std;

void updateX(VectorXd &theta, VectorXd &thetaPre, MatrixXd &X, unordered_set<int> &Alpha)
{
  double threadHold = 1e-8;
  double threadHoldc = 1e-5;
  int size = theta.size();
  int cur_rows_X = X.rows();
  int ncolX = X.cols();
  
  int  col = floor(sqrt(size));
  int  row = size/col;
  row = row +4;

  VectorXd Xi(ncolX);
  VectorXd Xi_tmp(ncolX);
  MatrixXd X_tmp = X; 
  

	// checking for new zeros in theta and updates X matrix accordingly
  // we may want to switch to use Alpha later on since working on theta
  // could potentially create duplicates. 
  for(int i = 0; i < size; i++)
  {
    if(theta(i) < threadHold &&  Alpha.find(i)== Alpha.end())
    {		
      Alpha.insert(i);
      Xi.fill(0);
			Xi(i) = -1;
			Xi_tmp(i) = Xi(i);
      X.conservativeResize(cur_rows_X+1,ncolX);
      X.row(cur_rows_X) = Xi;
      cur_rows_X +=1;
    }
  }
 

  for(int i = 1; i < col; i++)
  {
    for(int j = 1; j < col; j++)
    {
      if(abs(theta(i*col+j)-threadHoldc) < threadHold && Alpha.find(i*col+j)== Alpha.end() )
      {
		    Alpha.insert(i*col+j);
      	Xi.fill(0);
        Xi(i*col+j) = -1;
				Xi_tmp(i*col+j) = Xi(i*col+j) - 1;
        X.conservativeResize(cur_rows_X+1,ncolX);
        X.row(cur_rows_X) = Xi;
        cur_rows_X +=1;
	   }
    }
  }

	 
  // check for the sum of the all the I-spline equal to one or not
  if(abs(1-theta.sum()) <threadHold && Alpha.find(size)==Alpha.end())
  {
    Alpha.insert(size);
    Xi.fill(1);
    X.conservativeResize(cur_rows_X+1,ncolX);
    X.row(cur_rows_X) = Xi;
    cur_rows_X += 1;
  }

  return;	
}






