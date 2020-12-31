#include "RcppArmadillo.h"
// #include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

struct asset_info {
  double sum, sum2, stdev;
};

//[correlation matrix](http://en.wikipedia.org/wiki/Correlation_and_dependence).
// n,sX,sY,sXY,sX2,sY2
// cor = ( n * sXY - sX * sY ) / ( sqrt(n * sX2 - sX^2) * sqrt(n * sY2 - sY^2) )
inline asset_info compute_asset_info(const NumericMatrix& mat, 
				     const int icol, const int rstart, const int rend) {
  double sum, sum2;
  sum = sum2 = 0;

  for (int r = rstart; r < rend; r++) {
    double d = mat(r, icol);
    sum += d;
    sum2 += pow(d,2);
  }

  asset_info res;
  res.sum = sum;
  res.sum2 = sum2;
  res.stdev = sqrt((rend-rstart) * sum2 - pow(sum, 2));
  return res;
}



inline NumericMatrix c_cor_helper(const NumericMatrix& mat, const int rstart, const int rend) {
  int nc = mat.ncol();
  int nperiod = rend - rstart;
  NumericMatrix rmat(nc, nc);

  vector<asset_info> info(nc);
  for (int c = 0; c < nc; c++)
    info[c] = compute_asset_info(mat, c, rstart, rend);

  for (int c1 = 0; c1 < nc; c1++) {
    for (int c2 = 0; c2 < c1; c2++) {
      double sXY = 0;

      for (int r = rstart; r < rend; r++)
	sXY += mat(r, c1) * mat(r, c2);

      rmat(c1, c2) = (nperiod * sXY - info[c1].sum * info[c2].sum) / (info[c1].stdev * info[c2].stdev);
    }
  }

  return rmat;
}


/*

//' Fast selection of variables below correlation threshold
//' 
//' I'll need to fill this out.
//'
//' @param mat A matrix 
//' @param threshold The threshold. Should be a number between 0 and 1. 
//' @return A logical vector where the elements that are TRUE correspond to the variables that are not correlated above the threshold with any of the subsequent variables.
//' @author Claus Ekstrøm <claus@@rprimer.dk>
//' @examples
//'
//' x <- sample(10, 20, replace = TRUE)
//' bin(x, 15)
//' 
//' @export
// [[Rcpp::export]]
LogicalVector select_variables_using_cor(const NumericMatrix& mat, double threshold) {
  int rstart=0;
  int rend = mat.nrow();
  int nc = mat.ncol();
  int nperiod = rend - rstart;
  NumericMatrix rmat(nc, nc);
  LogicalVector keep(nc);

  vector<asset_info> info(nc);

  /*
  for (int c = 0; c < nc; c++)
    info[c] = compute_asset_info(mat, c, rstart, rend);

  for (int c1 = 0; c1 < nc; c1++) {
    for (int c2 = 0; c2 < c1; c2++) {
      double sXY = 0;

      for (int r = rstart; r < rend; r++)
	sXY += mat(r, c1) * mat(r, c2);

      rmat(c1, c2) = (nperiod * sXY - info[c1].sum * info[c2].sum) / (info[c1].stdev * info[c2].stdev);
    }
  }
  * /


  // From here
  for (int c = 0; c < nc; c++) {
    info[c] = compute_asset_info(mat, c, rstart, rend);
  }

  for (int c1 = 0; c1 < nc; c1++) {
    for (int c2 = 0; c2 < c1; c2++) {

      if (keep(c2))
        continue;

      double sXY = 0;

      // 
      for (int r = rstart; r < rend; r++)
	sXY += mat(r, c1) * mat(r, c2);

      

      if (fabs(nperiod * sXY - info[c1].sum * info[c2].sum) / (info[c1].stdev * info[c2].stdev) > threshold) {
        keep(c2)=1;
        break;
      }
    }
  }



  return (!keep);
}





//' Fast selection of variables below correlation threshold
//' 
//' I'll need to fill this out.
//'
//' @param mat A matrix 
//' @param threshold The threshold. Should be a number between 0 and 1. 
//' @return A logical vector where the elements that are TRUE correspond to the variables that are not correlated above the threshold with any of the subsequent variables.
//' @author Claus Ekstrøm <claus@@rprimer.dk>
//' @source https://arxiv.org/pdf/1810.11332.pdf
//' @examples
//'
//' x <- sample(10, 20, replace = TRUE)
//' bin(x, 15)
//' 
//' @export
// [[Rcpp::export]]
double dcov(const arma::vec & x, const arma::vec & y) {

  int n = x.size();

  // Sanity check

  arma::uvec indices = arma::sort_index(x);

  x = x(indices);
  y = y(indices);


  arma::vec si = cumsum(x);
}


*/

