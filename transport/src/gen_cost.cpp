#include <RcppEigen.h>
#ifdef _OPENMP
# include <omp.h>
#endif
using namespace Rcpp;
using namespace std;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::export]]

NumericMatrix gen_cost(NumericMatrix AR, NumericMatrix BR, int threads) {
  Eigen::Map<Eigen::MatrixXd> A(as<Eigen::Map<Eigen::MatrixXd> >(BR));
  Eigen::Map<Eigen::MatrixXd> B(as<Eigen::Map<Eigen::MatrixXd> >(AR));
  Eigen::setNbThreads(threads);
  int Al=A.rows();
  int Bl=B.rows();
  
  Eigen::MatrixXd x=(A*A.transpose()).diagonal();
  Eigen::MatrixXd z=(B*B.transpose()).diagonal();
  Eigen::MatrixXd onesA = Eigen::MatrixXd::Constant(1,Al, 1.0);
  Eigen::MatrixXd onesB = Eigen::MatrixXd::Constant(Bl,1, 1.0);
  Eigen::MatrixXd Cmat=(z*onesA+onesB*x.transpose()-(2.0*B*A.transpose()));
  
  
  return  wrap(Cmat);
}
