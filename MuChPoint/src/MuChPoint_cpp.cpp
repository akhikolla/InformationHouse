#include <Rcpp.h>
#include <cmath>
#include <iostream>

using namespace Rcpp;
using namespace std;

//' Compute the Delta of the dynamic programming
//'
//' Compute the Delta of the dynamic programming in \code{\link{Rcpp}}
//'
//' @param x the matrix of rank
//' @export Compute_Cn1n2
//' @rdname Compute_Cn1n2
// [[Rcpp::export]]
NumericMatrix Compute_Cn1n2(NumericMatrix x) {

  int n = x.nrow();
  NumericMatrix mat(x.nrow(),x.ncol());
  NumericVector vect(x.ncol());
  NumericVector vect_bis(x.ncol());
  NumericVector vect_suite(x.ncol());

  // calcul du vecteur des sommes de colonnes dans le cas de la diagonale
  for(int ii=0;ii<n;ii++)
  { double sum =0;
    for(int jj=0;jj<n;jj++){
      sum = sum + std::pow(x(jj,ii)/(2*n+2),2);
      vect(ii) = sum;
    }}



  for(int n1 = 0; n1 < n; ++n1) {
    for(int n2 = 0; n2 < n; ++n2) {

      // cas de la diagonale
      if(n1==n2) { mat(n1,n2) = vect(n1); }

      // cas courant
      else if (n2>n1) {

        double res = 0;
        for(int i=0;i< n;i++)
        {
          double sum2 = 0;
          for(int j=n1;j<=n2;j++){
            sum2 = sum2 + x(i,j);
          }
          vect_bis(i) = sum2;
          vect_suite(i) = std::pow(vect_bis(i)/((n2+1)-n1)-(n+1)/2.0,2);
          res=res + vect_suite(i)*((n2+1)-n1);
        }
        mat(n1,n2) = res;

      }

      // reste de la matrice remplie de zero
      else {
        mat(n1,n2)=0;
      }
      Rcpp::checkUserInterrupt();
    }
  }
  return mat;
}
