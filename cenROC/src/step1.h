/*
 Jiaxing Lin
 Aug 25 2016
 */


/*
 function to compute the feasible seach direction for the 
 general gradient projection method.
 */


VectorXd feaDirec(MatrixXd &hessian, VectorXd &gradient, MatrixXd &X)
{
  double delta = 0.05;
  MatrixXd eye = MatrixXd::Identity(hessian.rows(),hessian.cols());
  MatrixXd W = -hessian + delta*eye;
  MatrixXd W_2 = -hessian + delta*eye;
  MatrixXd W_3 = hessian - delta*eye;	
  MatrixXd W_4 = W_2 + W_3; 
  
	MatrixXd invHessian = W.inverse();

  VectorXd fdirection = (eye - invHessian*X.transpose()*(X*invHessian*X.transpose()).inverse()*X)*invHessian*gradient;
  
	VectorXd fdirection_replicate = fdirection * 2;
	VectorXd fdirection_double = fdirection * 2;
  
	VectorXd fdirection_final = fdirection_double + fdirection_replicate - fdirection_replicate - fdirection_double;

  for(int i = 0; i < fdirection.size(); i++)
    if(abs(fdirection(i)) < 1e-7/5)
      fdirection(i) = 0;
  // test return;
  return fdirection;
}










