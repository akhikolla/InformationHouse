
#ifndef _UTILS_CC_
#define _UTILS_CC_

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(.postmdiag)]]
NumericMatrix postmdiag(NumericMatrix X,NumericVector d)
{
  int c = X.ncol();
  int r = X.nrow();
  double dj;
  NumericMatrix Y(r,c);

  if(d.size() != c)
  {
    Rcpp::stop("Length of d must be same as number of columns of X");
  }
  for(int j=0;j<c;++j)
  {
    dj = d[j];
    NumericMatrix::iterator iX = X.begin()  + j*r;
    NumericMatrix::iterator iY = Y.begin()  + j*r;
    for(int i=0;i<r;++i)
      *(iY++) = *(iX++)*dj;
  }
  return(Y);
}


/*
 *The followuing function computes the 1/rowSum( (D*(X- 1*mu'))^2 ) for a dgCMatrix X
 * without additional storage
 */
// [[Rcpp::export(.icolSumSqdgC)]]
NumericVector colSumSqdgC(S4 X,NumericVector D,NumericVector mu)
{
  IntegerVector dims = X.slot("Dim");
  IntegerVector rowp = as<Rcpp::IntegerVector>(X.slot("i"));
  IntegerVector colp = as<Rcpp::IntegerVector>(X.slot("p"));
  NumericVector x = as<Rcpp::NumericVector>(X.slot("x"));

  int nrow = dims[0], ncol = dims[1];

  NumericVector v(ncol);
  double temp;
  for(int j=0;j<ncol;++j)
  {
    double vsum = 0.0;
    double muj = mu[j];
    int start = colp[j], end = colp[j+1];
    int iold = 0;
    for(int i=start; i<end;++i)
    {
      int irow = rowp[i];
      for(int l=iold;l<irow;++l)
      {
        temp = D[l]*muj;
        vsum += temp*temp;
      }
      temp = D[irow]*(x[i] - muj);
      vsum += temp*temp;
      iold = irow + 1;
    }
    for(int l=iold; l<nrow; ++l)
    {
      temp = D[l]*muj;
      vsum += temp*temp;
    }
    v[j] = 1.0/std::sqrt(vsum);
  }

  return(v);
}

/*
 * The following function computes the exact same thing as above but for regular matrices.
 *
 */
// [[Rcpp::export(.icolSumSq)]]
NumericVector colSumSq(NumericMatrix X,NumericVector D,NumericVector mu)
{
  int n = X.nrow(), p = X.ncol();
  NumericVector v(p);
  for(int j=0;j<p;++j)
  {
    double muj = mu[j];
    double vsum = 0.0;
    NumericMatrix::iterator ix = X.begin() + j*n;
    for(int i=0;i<n;++i)
    {
      double temp = D[i] * (*(ix++) - muj);
      vsum += temp*temp;
    }
    v[j] = 1/std::sqrt(vsum);
  }

  return(v);
}










#endif

