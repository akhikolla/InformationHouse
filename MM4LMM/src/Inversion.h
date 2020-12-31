#ifndef INVERSION
#define INVERSION
#include <RcppEigen.h>
#include <iostream>

using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;
void blocki(Eigen::MatrixXd & x, int x1, int n, Eigen::MatrixXd & y, int y1, double & log_det, double & det, double eps);

inline void sym_inverse(Eigen::MatrixXd & X,Eigen::MatrixXd & Y, double & log_det, double & det, double eps) {
  if (X.isDiagonal()) {
	VectorXd xDiag = X.diagonal();
	VectorXd xinv = xDiag.cwiseInverse();
	Y = xinv.asDiagonal();
	VectorXd logDiag(xDiag.size());
	for (int i = 0 ; i < xDiag.size(); i++) logDiag[i] = log(xDiag[i]);
	log_det = logDiag.sum();
	det = xDiag.prod();
  }else{
  	blocki(X, 0, X.rows(), Y, 0, log_det, det, eps);
  	Y.triangularView<Lower>() = Y.transpose(); // symétriser
  }
}

void blocki(Eigen::MatrixXf & x, int x1, int n, Eigen::MatrixXf & y, int y1, float & log_det, float & det, float eps);

inline void sym_inverse(Eigen::MatrixXf & X,Eigen::MatrixXf & Y, float & log_det, float & det, float eps) {
  blocki(X, 0, X.rows(), Y, 0, log_det, det, eps);
  Y.triangularView<Lower>() = Y.transpose(); // symétriser
}

template<typename T1, typename T2>
double trace_of_product(const Eigen::MatrixBase<T1> & a, const Eigen::MatrixBase<T2> & b) {
  int n(a.rows());
  int m(a.cols());
  if(b.rows() != m || a.cols() != n) {
    Rcerr << "\nSize mismatch in trace of product (" << n << "," << m << ") (" << b.rows() << "," << a.cols() << ")\n" ;
    return 0;
  }

  // sum i,k a_ik b_ki
  double S = 0;
  for(int i = 0; i < n; i++) {
    for(int k = 0; k < m; k++) {
      S += a(i,k)*b(k,i);
    }
  }
  return S;
}

#endif
