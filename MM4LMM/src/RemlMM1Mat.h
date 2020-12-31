// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
// #define ANY(_X_) (std::any_of(_X_.begin(), _X_.end(), [](bool x) {return x;})) 




#include <Rcpp.h>
#include "Inversion.h"

using namespace Rcpp;
using namespace Eigen;
typedef Map<MatrixXd> Map_MatrixXd;


List MM_Reml1MatRcpp(VectorXd & Y, MatrixXd & X, MatrixXd & VarInv , double logdetVar){

	int NbObs(Y.size());
	int NbCof(X.cols());

	double logdetxVarx , detxVarx;
	MatrixXd Varx = VarInv * X;
	MatrixXd xVarx = X.transpose() * Varx;
	MatrixXd xVarx_inv(NbCof,NbCof);
	sym_inverse(xVarx , xVarx_inv , logdetxVarx , detxVarx , 0);

	MatrixXd P = VarInv - Varx * xVarx_inv.selfadjointView<Lower>() * Varx.transpose() ;
  VectorXd InvFixed = P.selfadjointView<Lower>() * Y;

	
	double yPy = Y.transpose() * InvFixed;
	double SigmaNew =  yPy/(NbObs - NbCof);

	double LogLik = -(logdetxVarx + logdetVar + yPy)/2;
	
	VectorXd Beta = xVarx_inv.selfadjointView<Lower>() * (Varx.transpose() * Y);

	return List::create(Rcpp::Named("Beta")=Beta,Rcpp::Named("Sigma2")=SigmaNew,Rcpp::Named("VarBeta") = xVarx_inv*SigmaNew,Rcpp::Named("LogLik (Reml)")=LogLik,Rcpp::Named("NbIt")=1 , Rcpp::Named("Method")="Reml");
	
}                         
                          


