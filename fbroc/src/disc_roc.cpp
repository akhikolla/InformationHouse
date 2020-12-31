#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List add_roc_points(NumericVector tpr, NumericVector fpr) {
  int n = tpr.size();
  int segments = 1;
  IntegerVector in_segment (n);
  in_segment[0] = 1;
  // count segments and see in which segment each part is
  for (int i = 1; i < n; i++) {
    if ((tpr[i] != tpr[i-1]) && (fpr[i] != fpr[i-1])) segments++;
    in_segment[i] = segments;
  }
  
  List out(3);
  // if just one segment, basically do nothing
  if (segments == 1) {
    out[0] = tpr;
    out[1] = fpr;
    out[2] = in_segment;
    return out;
  }
  NumericVector tpr_new (n + (segments - 1));
  NumericVector fpr_new (n + (segments - 1));
  IntegerVector in_segment_new (n + (segments - 1));
  
  int current_index_new = 0;
  int current_index_old = 0;
  int current_segment = 1;
  
  while (current_index_new < (n + (segments - 1))) {
    while ((in_segment[current_index_old] == current_segment) && 
           (current_index_new < (n + (segments - 1)))) {
      tpr_new[current_index_new] = tpr[current_index_old];
      fpr_new[current_index_new] = fpr[current_index_old];
      in_segment_new[current_index_new] = in_segment[current_index_old];
      current_index_new++;
      current_index_old++;
    }
    if (current_index_new < (n + (segments - 1))) {
      tpr_new[current_index_new] = tpr[current_index_old];
      fpr_new[current_index_new] = fpr[current_index_old - 1];
      in_segment_new[current_index_new] = current_segment + 1;
      current_index_new++;
      current_segment = in_segment[current_index_old];
    }
  }
  
  out[0] = tpr_new;
  out[1] = fpr_new;
  out[2] = in_segment_new;
  return out;
}
