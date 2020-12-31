#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sumxxt(List x, int L) {
  List xlist(x);
  int n = xlist.size();

  std::vector<double> res((L * (L+1)) / 2);


  for(int l = 0; l < n; l++) {
    // Rcout<<"at line"<<l<<std::endl;

    SEXP ll = xlist[l];
    Rcpp::NumericVector xi(ll);

    for(int j = 0; j < L; j++)
      for(int i = j; i < L; i++) {
        // Rcout<<"adding "<<xi[i]<<"*"<<xi[j]<<"="<<xi[i] * xi[j]<<" to "<<(i * (i+1))/2 + j<<" so ("<<i<<", "<<j<<")"<<std::endl;
        res[(i * (i+1))/2 + j] += xi[i] * xi[j];
      }
  }

  return wrap(res);
}


/* this function takes a left position, right position, and elp, meant to be applied for each cluster
* it returns a vector which, for all unique time points, contains the sum of elp gathered from the cluster.
* i.e. say we have (0, 3) and (2,3), two individuals at risk for a total of 5 events (positions 0 to 4).
* Then this will give a vector with
* position 0 elp1
* position 1 elp1
* position 2 elp1 + elp2
* position 3 elp1 + elp2
* position 4 0.0
*/

// [[Rcpp::export]]
NumericVector cumsum_elp(NumericVector left, NumericVector right, NumericVector elp, int maxlength) {

  NumericVector x(maxlength, 0.0);
  unsigned int nrow = left.size();

  unsigned int upto = maxlength; // this could be the length of the out vector

  for(unsigned int i = 1; i <= upto; i++) {
    for(unsigned int j = 0; j < nrow; j++) {
      if((left[j] < i) && (i <= right[j]))
        x[i-1] += elp[j];
    }
  }

  return x;

}

// this one is faster I think but something is still not right with the margins of the looping
// //[[Rcpp::export]]
// NumericVector inf_mat_match(NumericVector left, NumericVector right, NumericVector summand, int ntimes) {
//
//   NumericVector x(ntimes, 0.0);
//
//   unsigned int nrow_max = left.size();
//   for(unsigned int nrow = 0; nrow < nrow_max; nrow++) {
//     // Rcout<<"at row..."<<nrow<<std::endl;
//     // Rcout<<"edges: left "<<left[nrow]<<" and right "<<right[nrow]<<std::endl;
//     for(int tp = left[nrow]; tp <= right[nrow]; tp++) {
//       // Rcout<<"x["<<tp<<"] += summand["<<nrow<<"]"<<std::endl;
//       x[tp] += summand[nrow];
//     }
//   }
//   return x;
//
// }

/*** R
# This function is used to calculate if you have a list of vectors of the same length
# and you want the sum of all the x%*%x transposed
# n <- 1
# L <- 3
# set.seed(1)
# x <- matrix(rnorm(n * L), nrow = n, ncol = L)
# xx <- apply(x, 1, list)
# xxx <- lapply(xx, function(x) x[[1]])
#
# xxx
#
# m <- matrix(0, L, L)
# a <- sumxxt(xxx, L)
# m[upper.tri(m, diag = TRUE)] <- a
#
# m3 <- Reduce("+", lapply(xxx, function(vec) vec %*% t(vec)))
#
# m2 <- m + t(m) - diag(diag(m))
# all.equal(m2, m3)
# all.equal(m[upper.tri(m)] , m3[upper.tri(m3)])


*/


// [[Rcpp::export]]
NumericVector rowsum_vec(NumericVector x, NumericVector pos, int lgth) {

  NumericVector res(lgth);

  for(int i = 0; i < x.size(); i++) {
    res[pos[i] - 1] += x[i];
  }
  for(int i = res.size() - 1; i!=0; i--) {
    res[i - 1] += res[i];
  }

  return(res);
}


/*** R
# take it from the cgd
# data(cgd)
#
#   ord1 <- match(cgd$tstop, sort(unique(cgd$tstop)))
#   library(microbenchmark)
#
#   all.equal(
#     with(cgd, rowsum_vec(elp, ord1, max(ord1))),
#     with(cgd, rev(cumsum(rev(rowsum(elp, tstop)))))
#   )
#
#   microbenchmark(
#     with(cgd, rowsum_vec(elp, ord1, max(ord1))),
#     with(cgd, rev(cumsum(rev(rowsum(elp, tstop)))))
#   )
#
# #with the rats
# # take it from the cgd
#   data(rats)
#   head(rats)
#   rats$elp <- rnorm(nrow(rats))
#   ord1 <- match(rats$time, sort(unique(rats$time)))
#   library(microbenchmark)
#
#   all.equal(
#     with(rats, rowsum_vec(elp, ord1, max(ord1))),
#     with(rats, rev(cumsum(rev(rowsum(elp, time)))))
#   )
#
#   microbenchmark(
#     with(rats, rowsum_vec(elp, ord1, max(ord1))),
#     with(rats, rev(cumsum(rev(rowsum(elp, time)))))
#   )

# 10x reduction in running time

  */

