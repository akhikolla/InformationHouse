#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

 using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix k_chol(NumericMatrix GlobalMat)
 {
  NumericMatrix Mat = clone(GlobalMat);

  char uplo = 'U';
  int n = Mat.nrow();
  int lda = n;
  int info = 0;
  for (int i = 0; i < n; i++) 	/* zero the lower triangle */
	for (int j = i+1; j < n; j++) Mat[j + n * i] = 0.;

  F77_CALL(dpotrf)(&uplo, &n, Mat.begin(), &lda, &info);
  return Mat;
}

//[[Rcpp::export]]
NumericMatrix k_chol2inv(NumericMatrix GlobalMat)
 {
  NumericMatrix Mat = clone(GlobalMat);

  char uplo = 'U';
  int n = Mat.nrow();
  int lda = n;  //std::max(1, n)
  int info = 0;

  F77_CALL(dpotri)(&uplo, &n, Mat.begin(), &lda, &info);

  for (int i = 0; i < n; i++)
      for (int j = i+1; j < n; j++)
    Mat[j + i * n] = Mat[i + j * n];

  return Mat;
}

//[[Rcpp::export]]
NumericVector k_eigenvalue(NumericMatrix GlobalMat)
{
  NumericMatrix Mat = clone(GlobalMat); //Work on a copy
  char jobz = 'N', range = 'A',  uplo = 'U';
  int n = Mat.nrow(), info = 0, lda = n, iu = n;
  int il = 1; //not referenced if range = 'A'
  double vl = 0.0, vu = 0.0;
  double abstol = 1.490116e-08; // tolerance: 0.0 or 1.490116e-08
  int m = n, ldz = n;
  IntegerVector isuppz(2*m);
  int lwork = -1, liwork = -1;
  double work_tmp;
  int iwork_tmp;
  NumericVector lambda(n);
  NumericMatrix eigvec(ldz, m);

  // Query for the optimal work array.
  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, Mat.begin(), &lda, &vl, &vu, &il, &iu, &abstol,
           &m, lambda.begin(), eigvec.begin(), &ldz, isuppz.begin(), &work_tmp, &lwork, &iwork_tmp, &liwork, &info);

  lwork = work_tmp;
  liwork = iwork_tmp;
  NumericVector work(lwork);
  IntegerVector iwork(liwork);

  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, Mat.begin(), &lda, &vl, &vu, &il, &iu, &abstol,
           &m, lambda.begin(), eigvec.begin(), &ldz, isuppz.begin(), work.begin(), &lwork, iwork.begin(), &liwork, &info);

  std::reverse(lambda.begin(), lambda.end()); // Reverse the order
  return lambda;
}

//[[Rcpp::export]]
NumericMatrix k_inv(NumericMatrix GlobalMat)
 {
  NumericMatrix Mat = clone(GlobalMat);

  char uplo = 'U';
  int n = Mat.nrow();
  int lda = n;  //std::max(1, n)
  int info = 0;

  F77_CALL(dpotrf)(&uplo, &n, Mat.begin(), &lda, &info);
  F77_CALL(dpotri)(&uplo, &n, Mat.begin(), &lda, &info);

  for (int i = 0; i < n; i++)
      for (int j = i+1; j < n; j++)
    Mat[j + i * n] = Mat[i + j * n];

  return Mat;
}

/* Euclidean distance function*/
double dist(double east1, double north1, double east2, double north2) {
  double e_distance = sqrt(pow(east1 - east2, 2) + pow(north1 - north2, 2)) ;
  return e_distance;
}

//[[Rcpp::export]]
NumericMatrix k_distmat (NumericMatrix Mat)
  {
  int N = Mat.nrow();
  NumericMatrix DistMat(N, N);
  for (int i = 0 ; i < N ; i++) {
    for (int j = 0 ; j < N ; j ++) {
      DistMat(i, j) = dist(Mat(i, 0),Mat(i, 1),Mat(j, 0),Mat(j, 1));
    }
  }
  return DistMat;
}
