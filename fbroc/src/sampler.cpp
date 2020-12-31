#include <Rcpp.h>
using namespace Rcpp;
#include "sampler.h"


IntegerVector Sampler_base::get_shuffled_index(bool which_class) const
{
  if (which_class) return shuffled_pos_index;
  else return shuffled_neg_index;
}

Sampler_base::~Sampler_base() {
  // do nothing
}

Sampler_Stratified::~Sampler_Stratified() {
  // do nothing
}

Sampler_Stratified::Sampler_Stratified(IntegerVector true_class)
{
  n = true_class.size();

  n_pos = 0;
  n_neg = 0;
  for (int i = 0; i < n; i++) {
    if (true_class[i] == 1) {
      n_pos++;
    } else n_neg++;
  }
  
  shuffled_pos_index = IntegerVector(n_pos);
  shuffled_neg_index = IntegerVector(n_neg);
}

void Sampler_Stratified::generate()
{
  // call R for random number generation
  NumericVector random_pos = runif(n_pos);
  NumericVector random_neg = runif(n_neg);
  
  // use stratified bootstrapping
  for (int i = 0; i < n_pos; i++) {
    shuffled_pos_index[i] = (int)(n_pos * random_pos[i]);
  }
  for (int i = 0; i < n_neg; i++) {

    shuffled_neg_index[i] = (int)(n_neg * random_neg[i]);
  }
}
