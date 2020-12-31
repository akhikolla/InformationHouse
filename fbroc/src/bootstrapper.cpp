#include <Rcpp.h>
using namespace Rcpp;

#include "roc.h"
#include "sampler.h"
#include "bootstrapper.h"


void Bootstrapped_ROC::bootstrap() {
  sampler->generate();
  IntegerVector shuffled_pos_index = sampler->get_shuffled_index(true);
  IntegerVector shuffled_neg_index = sampler->get_shuffled_index(false);
  strat_shuffle(shuffled_pos_index, shuffled_neg_index);
}

Bootstrapped_ROC::Bootstrapped_ROC(NumericVector pred, IntegerVector true_class):ROC(pred, true_class)
{
  sampler = new Sampler_Stratified(true_class);
}

Bootstrapped_ROC::~Bootstrapped_ROC()
{
  delete sampler;
}
