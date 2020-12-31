#include <Rcpp.h>
using namespace Rcpp;

#include "roc.h"
#include "sampler.h"
#include "paired.roc.h"


void Bootstrapped_paired_ROC::bootstrap() {
  sampler->generate();
  IntegerVector shuffled_pos_index = sampler->get_shuffled_index(true);
  IntegerVector shuffled_neg_index = sampler->get_shuffled_index(false);
  roc[0].strat_shuffle(shuffled_pos_index, shuffled_neg_index);
  roc[1].strat_shuffle(shuffled_pos_index, shuffled_neg_index);
}


ROC & Bootstrapped_paired_ROC::get_roc(int number) {
  return roc[number];
};

Bootstrapped_paired_ROC::Bootstrapped_paired_ROC(NumericVector pred1, 
                                                 NumericVector pred2,
                                                 IntegerVector true_class)
{ 
  roc[0] = ROC(pred1, true_class);
  roc[1] = ROC(pred2, true_class);
  sampler = new Sampler_Stratified(true_class);
}

Bootstrapped_paired_ROC::~Bootstrapped_paired_ROC()
{
  delete sampler;
}
