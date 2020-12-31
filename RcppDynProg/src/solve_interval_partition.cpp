
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::IntegerVector;



IntegerVector solve_interval_partition_k_worker(const NumericMatrix &x, int kmax) {
  // for cleaner notation
  // solution and x will be indexed from 1 using
  // R_INDEX_DELTA
  // intermediate arrays will be padded so indexing
  // does not need to be shifted
  const int R_INDEX_DELTA = -1;
  const int R_SIZE_PAD = 1;
  
  // get shape of problem
  const int n = x.nrow();
  if(kmax>n) {
    kmax = n;
  }
  
  // get some edge-cases
  if((kmax<=1)||(n<=1)) {
    IntegerVector solution(2);
    solution(1 + R_INDEX_DELTA) = 1;
    solution(2 + R_INDEX_DELTA) = n+1;
    return solution;
  }
  
  if(n != x.ncol()) {
    throw std::range_error("Inadmissible value");
  }
  
  // best path cost up to i (row) with exactly k-steps (column)
  arma::Mat<double> path_costs(n + R_SIZE_PAD, 2, arma::fill::none);
  // how many steps we actually took
  arma::Mat<int> k_actual(n + R_SIZE_PAD, 2, arma::fill::ones);
  // how we realized each above cost
  arma::Mat<int> prev_step(n + R_SIZE_PAD, kmax + R_SIZE_PAD, arma::fill::ones);
  
  // fill in path and costs tables
  for(int i=1; i<=n; ++i) {
    const double xi = x(1 + R_INDEX_DELTA, i + R_INDEX_DELTA);
    path_costs(i, 0) = xi;
    path_costs(i, 1) = xi;
  }
  // refine dynprog table
  int kcurrent = 0;
  for(int ksteps=2; ksteps<=kmax; ++ksteps) {
    kcurrent = ksteps%2;
    const int kprev = 1 - kcurrent;
    // compute larger paths
    for(int i=1; i<=n; ++i) {
      // no split case
      int pick = i;
      int k_seen = 1;
      double pick_cost = x(1 + R_INDEX_DELTA, i + R_INDEX_DELTA);
      // split cases
      for(int candidate=1; candidate<i; ++candidate) {
        const double cost = path_costs(candidate, kprev) +
          x(candidate + 1 + R_INDEX_DELTA, i + R_INDEX_DELTA);
        const int k_cost = k_actual(candidate, kprev) + 1;
        if((cost<=pick_cost) &&
           ((cost<pick_cost)||(k_cost<k_seen))) {
          pick = candidate;
          pick_cost = cost;
          k_seen = k_cost;
        }
      }
      path_costs(i, kcurrent) = pick_cost;
      prev_step(i, ksteps) = pick;
      k_actual(i, kcurrent) = k_seen;
    }
  }
  
  // now back-chain for solution
  const int k_opt = k_actual(n, kcurrent);
  IntegerVector solution(k_opt+1);
  solution(1 + R_INDEX_DELTA) = 1;
  solution(k_opt + 1 + R_INDEX_DELTA) = n+1;
  int i_at = n;
  int k_at = k_opt;
  while(k_at>1) {
    const int prev_i = prev_step(i_at, k_at);
    solution(k_at + R_INDEX_DELTA) = prev_i + 1;
    i_at = prev_i;
    k_at = k_at - 1;
  }
  
  return solution;
}



//' solve_interval_partition interval partition problem with a bound on number of steps.
//' 
//' Solve a for a minimal cost partition of the integers [1,...,nrow(x)] problem where for j>=i x(i,j).
//' is the cost of choosing the partition element [i,...,j]. 
//' Returned solution is an ordered vector v of length k<=kmax where: v[1]==1, v[k]==nrow(x)+1, and the 
//' partition is of the form [v[i], v[i+1]) (intervals open on the right).
//' 
//' @param x square NumericMatix, for j>=i x(i,j) is the cost of partition element [i,...,j] (inclusive).
//' @param kmax int, maximum number of segments in solution. 
//' @return dynamic program solution.
//' 
//' @examples
//' 
//' costs <- matrix(c(1.5, NA ,NA ,1 ,0 , NA, 5, -1, 1), nrow = 3)
//' solve_interval_partition(costs, nrow(costs))
//'
//' @export
// [[Rcpp::export]]
IntegerVector solve_interval_partition_k(NumericMatrix x, int kmax) {
  return solve_interval_partition_k_worker(x, kmax);
}








IntegerVector solve_interval_partition_no_k_worker(const NumericMatrix &x) {
  // for cleaner notation
  // solution and x will be indexed from 1 using
  // R_INDEX_DELTA
  // intermediate arrays will be padded so indexing
  // does not need to be shifted
  const int R_INDEX_DELTA = -1;
  const int R_SIZE_PAD = 1;
  
  // get shape of problem
  const int n = x.nrow();
  
  // get some edge-cases
  if(n<=1) {
    IntegerVector solution(2);
    solution(1 + R_INDEX_DELTA) = 1;
    solution(2 + R_INDEX_DELTA) = n+1;
    return solution;
  }
  
  if(n != x.ncol()) {
    throw std::range_error("Inadmissible value");
  }
  
  // best path cost up to i (row) with exactly k-steps (column)
  arma::Col<double> path_costs(n + R_SIZE_PAD, arma::fill::zeros);
  // how we realized each above cost entry is rhs of last interval (0 at start)
  arma::Col<int> prev_step(n + R_SIZE_PAD, arma::fill::ones);
  
  // fill in path and costs tables
  // refine dynprog table
  // compute larger paths
  for(int i=1; i<=n; ++i) {
    // no split case
    int pick = 0;
    double pick_cost = x(1 + R_INDEX_DELTA, i + R_INDEX_DELTA);
    // split cases
    for(int candidate=1; candidate<i; ++candidate) {
      const double cost = path_costs(candidate) +
        x(candidate + 1 + R_INDEX_DELTA, i + R_INDEX_DELTA);
      if(cost<=pick_cost) {
        pick = candidate;
        pick_cost = cost;
      }
    }
    path_costs(i) = pick_cost;
    prev_step(i) = pick;
  }
  
  // now back-chain for solution
  int k_at = 1;
  int i_at = n;
  while(prev_step(i_at)>0) {
    k_at = k_at + 1;
    i_at = prev_step(i_at);
  }
  const int k_opt = k_at;
  IntegerVector solution(k_opt+1);
  solution(1 + R_INDEX_DELTA) = 1;
  solution(k_opt + 1 + R_INDEX_DELTA) = n+1;
  i_at = n;
  k_at = k_opt;
  while(k_at>1) {
    const int prev_i = prev_step(i_at);
    solution(k_at + R_INDEX_DELTA) = prev_i + 1;
    i_at = prev_i;
    k_at = k_at - 1;
  }
  
  return solution;
}

//' solve_interval_partition interval partition problem, no boun on the number of steps.
//' 
//' Not working yet.
//' 
//' Solve a for a minimal cost partition of the integers [1,...,nrow(x)] problem where for j>=i x(i,j).
//' is the cost of choosing the partition element [i,...,j]. 
//' Returned solution is an ordered vector v of length k where: v[1]==1, v[k]==nrow(x)+1, and the 
//' partition is of the form [v[i], v[i+1]) (intervals open on the right).
//' 
//' @param x square NumericMatix, for j>=i x(i,j) is the cost of partition element [i,...,j] (inclusive).
//' @return dynamic program solution.
//' 
//' @examples
//' 
//' costs <- matrix(c(1.5, NA ,NA ,1 ,0 , NA, 5, -1, 1), nrow = 3)
//' solve_interval_partition(costs, nrow(costs))
//'
//' @export
// [[Rcpp::export]]
IntegerVector solve_interval_partition_no_k(NumericMatrix x) {
  return solve_interval_partition_no_k_worker(x);
}






//' solve_interval_partition interval partition problem.
//' 
//' Solve a for a minimal cost partition of the integers [1,...,nrow(x)] problem where for j>=i x(i,j).
//' is the cost of choosing the partition element [i,...,j]. 
//' Returned solution is an ordered vector v of length k<=kmax where: v[1]==1, v[k]==nrow(x)+1, and the 
//' partition is of the form [v[i], v[i+1]) (intervals open on the right).
//' 
//' @param x square NumericMatix, for j>=i x(i,j) is the cost of partition element [i,...,j] (inclusive).
//' @param kmax int, maximum number of segments in solution. 
//' @return dynamic program solution.
//' 
//' @examples
//' 
//' costs <- matrix(c(1.5, NA ,NA ,1 ,0 , NA, 5, -1, 1), nrow = 3)
//' solve_interval_partition(costs, nrow(costs))
//'
//' @export
// [[Rcpp::export]]
IntegerVector solve_interval_partition(NumericMatrix x, const int kmax) {
  IntegerVector soln1 = solve_interval_partition_no_k_worker(x);
  if(soln1.length()<=(kmax+1)) {
    return(soln1);
  }
  return(solve_interval_partition_k_worker(x, kmax));
}


