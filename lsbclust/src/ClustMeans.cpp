#include <Rcpp.h>
using namespace Rcpp;

//' @rdname ClustMeans
//' @title C++ Function for Cluster Means
//' @name ClustMeans
//' @description This function calculates the cluster means in vectorized form based on the current
//' value of the clustering vector.
//' @param nclust The number of clusters.
//' @param start The current clustering vector.
//' @param data The concatenated data, with J * K rows and N columns
//' @return A numeric matrix with \code{nclust} rows and \code{J*K} columns.
// [[Rcpp::export]]
NumericMatrix ClustMeans(int nclust, IntegerVector start, NumericMatrix data) {
   int n = data.ncol();
   int JK = data.nrow();
   
   int clust; 
   NumericVector nvec(nclust);
   NumericMatrix out(nclust, JK);

   for (int i = 0; i < n; i++) {
     clust = start(i) - 1;
     nvec(clust) += 1;
     out(clust,_) = out(clust,_) + data(_,i);
   }
   
   for (int i = 0; i < nclust; i++) {
     out(i,_) = out(i,_)/nvec(i);
   }

   return out;
//   return (Rcpp::List::create(Rcpp::Named("means")= out, Rcpp::Named("nvec") = nvec));

}
