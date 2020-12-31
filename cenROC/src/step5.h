/*
 Jiaxing Lin
 Aug 25 2016
 */

/*
 Checking for the stopping criterion
 */
#include <iostream>
#include <unordered_set>

void removeRow(MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();
  
  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols)
    = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
  
  matrix.conservativeResize(numRows,numCols);
}

VectorXd computeLambda(MatrixXd &hessian, VectorXd &gradient, MatrixXd &X)
{
  // compute the inversion of hessian matrix, which seem unevidable
  double delta = 0.05;
  MatrixXd eye = MatrixXd::Identity(hessian.rows(),hessian.cols());
  MatrixXd W = -hessian + delta*eye;

  MatrixXd invHessian = W.inverse();

  MatrixXd XInvHessXt = X*invHessian*X.transpose();
  MatrixXd invXInvHessXt = XInvHessXt.inverse();
  VectorXd lambda = invXInvHessXt*X*invHessian*gradient;
  return lambda;
}

void updateX_step5(VectorXd &lambda, MatrixXd &X, unordered_set<int> &Alpha)
{
  int jStar = 0;
  double minLamd = 1e100;
  for(int i = 0; i < lambda.size(); i++){
    if(lambda(i) < 0 && minLamd > lambda(i))
    {
      minLamd = lambda(i);
      jStar  	= i;
    }
    
  }
  
  int index = 0;
  if(X.row(jStar).sum()==X.cols())
    index = X.cols()+1;
  for(int i = 0; i < X.cols(); i++)
    if(X(jStar,i)== -1)
      index = i;
    removeRow(X,jStar);		
    Alpha.erase(index);
  return;
}

bool checkAllPos(VectorXd &lambda)
{
  int i = 0;
  while(i < lambda.size())
  {
    if(lambda(i) >= 0)
      i++;
    else
      return false;
  }
  return true;
}

bool checkConvergeAndUpdateX(VectorXd &fdirection, double epsilon, 
                             VectorXd &theta, MatrixXd &X, MatrixXd &hessian, 
                             VectorXd &gradient, unordered_set<int> &Alpha)
{
  double dnorm = fdirection.norm();
  if(dnorm > epsilon)
    return false;
  else{
    VectorXd lambda = computeLambda(hessian, gradient, X);
    if(checkAllPos(lambda))
      return true;	
     		else{
    			updateX_step5(lambda, X, Alpha);
     			return false;
     		}

  }
}



