#include <Rcpp.h>
using namespace Rcpp;

//' Convolution function for semi-Markov models
//' 
//' @description Internal function of 
//' \code{probtrans_by_convolution_semiMarkov}. It is written in C++
//' and is not meant to be called
//' directly by the user.
//' 
//' @param time_vector,integrand_1,integrand_2 Numeric vectors.
//' @author Moritz Gerstung & Rui Costa
//' @seealso \code{\link{probtrans_by_convolution_semiMarkov}}.
// [[Rcpp::export]]
NumericVector convolute_semiMarkov(NumericVector time_vector, NumericVector integrand_1, NumericVector integrand_2 ) {
  int xLen = time_vector.size();
  NumericVector overallSurvival(xLen);
  for(int i = 0; i < xLen; ++i) overallSurvival[i] = 0;
  for(int i = 1; i < xLen; ++i){
    for(int j = 1; j < i; ++j){
      overallSurvival[i] += integrand_2[i-j]*integrand_1[j];
    }
  }
  return overallSurvival;
}
//' Convolution function for Markov models
//' 
//' @description Internal function of 
//' \code{probtrans_by_convolution_Markov}. It is written in C++ and is
//'  not meant to be called
//' directly by the user.
//' 
//' @param time_vector,diff_vector,probtrans_vector_1,probtrans_vector_2 Numeric vectors.
//' @author Moritz Gerstung & Rui Costa
//' @seealso \code{\link{probtrans_by_convolution_Markov}}.
// [[Rcpp::export]]
NumericVector convolute_Markov(NumericVector time_vector, NumericVector diff_vector, NumericVector probtrans_vector_1, NumericVector probtrans_vector_2 ) {
  int xLen = time_vector.size();
  NumericVector overallSurvival(xLen);
  for(int i = 0; i < xLen; ++i) overallSurvival[i] = 0;
  for(int i = 1; i < xLen; ++i){
    for(int j = 1; j < i; ++j){
      overallSurvival[i] += (probtrans_vector_2[i]/probtrans_vector_2[j])*diff_vector[j] *probtrans_vector_1[j];
    }
  }
  return overallSurvival;
}
