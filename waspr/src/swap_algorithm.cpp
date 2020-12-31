// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

//' The swapping algorithm for computing Wasserstein barycenters
//'
//' @name swap_rcpp
//'
//' @param samples A cube containing samples for all subset posteriors (rows =
//'   subsets, columns = par, slices = samples)
//' @param acc accuracy
//' @param iter maximum number of iterations of the algorithm
//' @param out boolean indicating whether output for each iteration should be displayed (default = false)
//'
//' @return a three dimensional array (rows = subsets, columns = par, slices =
//'   samples) containing output from the swapping algorithm.
//'
// [[Rcpp::export]]

arma::cube swap_rcpp(arma::cube samples, double acc = 0.001, int iter = 10, bool out = false) {

  arma::cube samples_new = cube(samples);

  int d = samples.n_cols; //amount of parameters
  int n = samples.n_rows; //amount of subposteriors
  int k = samples.n_slices; //amount of posterior samples

  arma::cube cost(n, k, k);
  double cost_prev = -2*acc; //Assign previous cost such that stopping condition is never met the first iteration
  double cost_cur = 0;
  arma::vec cost_vec = arma::zeros(k);
  int its = 0;

  for (int i = 0; i < n; ++i){

   arma::cube samp_temp = samples_new;  //create temporary samples cube
   arma::mat samp = samples_new.row(i); //extract one subposterior
   samp_temp.shed_row(i);
   arma::mat sumsamples = sum(samp_temp, 0); //sum values of other subposteriors

    for (int k1 = 0; k1 < k; ++k1){
      for (int k2 = 0; k2 < k; ++k2){

      cost(i, k1, k2) = dot(samp.col(k1), sumsamples.col(k2));

      }
    }

    arma::mat cost_temp = cost.row(i);
    cost_vec += arma::diagvec(cost_temp);

  }

  cost_cur = mean(cost_vec);

  while ((std::abs(cost_cur - cost_prev) > acc) & (its < iter)){

    cost_prev = cost_cur; //assign current cost to previous
    cost_cur = 0; //set current cost to 0
    its += 1; //increase iteration counter by one
    cost_vec = arma::zeros(k);

    if(out) {
      Rcout << "Iteration:" << its << "\n";
      Rcout << "Cost:" << cost_prev << "\n";
    }

    for (int i = 0; i < n; ++i){

      arma::cube samp_temp = samples_new;  //create temporary samples cube
      arma::mat samp = samples_new.row(i); //extract one subposterior
      samp_temp.shed_row(i);
      arma::mat sumsamples = sum(samp_temp, 0); //sum values of other subposteriors

      for (int k1 = 0; k1 < k; ++k1){
        for (int k2 = 0; k2 < k; ++k2){

          //compute inner product for all combinations of posterior samples
          cost(i, k1, k2) = dot(samp.col(k1), sumsamples.col(k2));

          //compute inner product for specific combinations of posterior samples
          double inner_11 = dot(samp.col(k1), sumsamples.col(k1));
          double inner_12 = dot(samp.col(k1), sumsamples.col(k2));
          double inner_21 = dot(samp.col(k2), sumsamples.col(k1));
          double inner_22 = dot(samp.col(k2), sumsamples.col(k2));
          double left = inner_11 + inner_22;
          double right = inner_12 + inner_21;

          //Check swapping condition and swap
          if((right-left) > 0){
            arma::rowvec current_k1 = samples_new.subcube(i,0, k1, i, d-1, k1);
            samples_new.subcube(i,0, k1, i, d-1, k1) = samples_new.subcube(i,0, k2, i, d-1, k2);
            samples_new.subcube(i,0, k2, i, d-1, k2) = current_k1.t();
          }
        }
      }

      arma::mat cost_temp = cost.row(i);
      cost_vec += arma::diagvec(cost_temp); //compute current costs (sum over subposteriors)

    }

    cost_cur = mean(cost_vec); //compute current costs (mean over posterior samples)

  }



  return(samples_new); //Return rearranged subposteriors


}



