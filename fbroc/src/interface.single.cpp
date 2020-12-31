#include <Rcpp.h>
using namespace Rcpp;

#include "roc.h"
#include "sampler.h"
#include "bootstrapper.h"
#include "interface.common.h"
#include "performance.h"


// [[Rcpp::export]]
List tpr_fpr_boot(NumericVector pred, IntegerVector true_class, int n_boot) {

  Bootstrapped_ROC boot_roc (pred, true_class);
  int n_thres = boot_roc.get_n_thres();
  NumericMatrix tpr (n_boot, n_thres);
  NumericMatrix fpr (n_boot, n_thres);
  for (int i = 0; i < n_boot; i++) {
    boot_roc.bootstrap();
    tpr(i, _) = boot_roc.get_tpr();
    fpr(i, _) = boot_roc.get_fpr();
  }
  List out(2);
  out[0] = tpr;
  out[1] = fpr;
  return out;
}

// [[Rcpp::export]]
NumericVector get_uncached_perf(NumericVector pred, IntegerVector true_class, NumericVector param,
                                int n_boot, int measure) {
  PerfFun choosen_measure = pick_measure(static_cast <Measure>(measure));                                
  Bootstrapped_ROC boot_roc (pred, true_class);
  NumericVector roc_perf = NumericVector (n_boot);
  for (int i = 0; i < n_boot; i++) {
    boot_roc.bootstrap();
    roc_perf[i] = choosen_measure(boot_roc.get_tpr(), boot_roc.get_fpr(), param);
  }
  return roc_perf;
}

// [[Rcpp::export]]
NumericVector get_cached_perf(NumericMatrix tpr, NumericMatrix fpr, NumericVector param, int measure) {
  PerfFun choosen_measure = pick_measure(static_cast <Measure>(measure));
  int n_boot = tpr.nrow();
  NumericVector roc_perf = NumericVector (n_boot);
  //iterate over bootstrap replicates and get performance for each   
  for (int i = 0; i < n_boot; i++) {
    NumericVector tpr_vec = tpr(i, _);
    NumericVector fpr_vec = fpr(i, _);
    double perf = choosen_measure(tpr_vec, fpr_vec, param);
    roc_perf[i] = perf;
  }
  
  return roc_perf;
}

// [[Rcpp::export]]
NumericMatrix tpr_at_fpr_uncached(NumericVector pred, IntegerVector true_class, int n_boot, int n_steps) {
  Bootstrapped_ROC boot_roc (pred, true_class);
  NumericVector steps = get_steps(n_steps);
  NumericMatrix tpr_matrix (n_boot, n_steps + 1);
  for (int j = 0; j < n_boot; j++) { 
    boot_roc.bootstrap();
    tpr_matrix(j, _) = boot_roc.get_tpr_at_fpr(steps);
  }
  return tpr_matrix;
}

// [[Rcpp::export]]
NumericMatrix tpr_at_fpr_cached(NumericMatrix tpr, NumericMatrix fpr, int n_steps) {
  NumericVector steps = get_steps(n_steps);
  int n_boot = tpr.nrow();
  NumericMatrix tpr_matrix (n_boot, n_steps + 1);
  for (int j = 0; j < n_boot; j++) {
     NumericVector tpr_v = tpr(j, _);
     NumericVector fpr_v = fpr(j, _);
     tpr_matrix(j, _) = ROC::get_tpr_at_fpr(tpr_v, fpr_v, steps);
  }
  return tpr_matrix;
}

// [[Rcpp::export]]
NumericMatrix fpr_at_tpr_cached(NumericMatrix tpr, NumericMatrix fpr, int n_steps) {
  NumericVector steps = get_steps(n_steps);
  int n_boot = fpr.nrow();
  NumericMatrix fpr_matrix (n_boot, n_steps + 1);
  for (int j = 0; j < n_boot; j++) {
    NumericVector tpr_v = tpr(j, _);
    NumericVector fpr_v = fpr(j, _);
    fpr_matrix(j, _) = ROC::get_fpr_at_tpr(tpr_v, fpr_v, steps);
  }
  return fpr_matrix;
}

// [[Rcpp::export]]
NumericMatrix fpr_at_tpr_uncached(NumericVector pred, IntegerVector true_class, int n_boot, int n_steps) {
  Bootstrapped_ROC boot_roc (pred, true_class);
  NumericVector steps = get_steps(n_steps);
  NumericMatrix fpr_matrix (n_boot, n_steps + 1);
  for (int j = 0; j < n_boot; j++) { 
    boot_roc.bootstrap();
    fpr_matrix(j, _) = boot_roc.get_fpr_at_tpr(steps);
  }
  
  return fpr_matrix;
}

