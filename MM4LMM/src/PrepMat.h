// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
// #define ANY(_X_) (std::any_of(_X_.begin(), _X_.end(), [](bool x) {return x;})) 




#include <Rcpp.h>


using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;



List PrepMatRcpp(VectorXd Y , MatrixXd K1 , MatrixXd K2){
	GeneralizedSelfAdjointEigenSolver<MatrixXd> es(K1,K2);
	
	return List::create(Rcpp::Named("Diag")=es.eigenvalues(),Rcpp::Named("U")=es.eigenvectors(),Rcpp::Named("Ytilde")=es.eigenvectors().transpose() * Y);
}

