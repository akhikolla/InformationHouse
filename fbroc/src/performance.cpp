#include <Rcpp.h>
using namespace Rcpp;

#include "performance.h"
#include <algorithm>

double get_perf_auc(NumericVector &tpr, NumericVector &fpr, NumericVector &param) 
{
  int n_thres = tpr.size();
  double auc = 0.;
  // Numerical integration of step functions is easy
  for (int j = 1; j < n_thres; j++) {
    auc += (tpr[j - 1] - tpr[j]) * (2 - fpr[j - 1] - fpr[j]);
  }
  auc = 0.5 * auc;
  return auc;
}

double pauc_fpr_area(double fpr, NumericVector &param)
{
  if (fpr > param[1]) return 0;
  if (fpr < param[0]) return (param[1] - param[0]);
  return (param[1] - fpr);
}

double pauc_tpr_area(NumericVector &tpr, NumericVector &fpr, NumericVector &param, int index)
{
  if (tpr(index - 1) == tpr[index]) return 0; // necessary check to avoid division by zero later
  if (tpr[index - 1] < param[0]) return 0;
  if (tpr[index] > param[1]) return 0;
  
  double left = std::max(tpr[index], param[0]);
  double right = std::min(tpr[index - 1], param[1]);
  
  double base_val = 1 - fpr[index];
  double slope = (fpr[index] - fpr[index - 1]) / (tpr[index - 1] - tpr[index]);

  double value_left = base_val + (left - tpr[index]) * slope;
  double value_right = base_val + (right - tpr[index]) * slope;
  
  return (right - left) * (value_left + value_right);
}

double get_perf_pauc_over_fpr(NumericVector &tpr, NumericVector &fpr, NumericVector &param) 
{
  int n_thres = tpr.size();
  double p_auc = 0.;
  // Numerical integration of step functions is easy
  for (int j = 1; j < n_thres; j++) {
    p_auc += (tpr[j - 1] - tpr[j]) * (pauc_fpr_area(fpr[j - 1], param) 
                                      + pauc_fpr_area(fpr[j], param));
  }
  p_auc = 0.5 * p_auc;
  return p_auc;
}

double get_perf_pauc_over_tpr(NumericVector &tpr, NumericVector &fpr, NumericVector &param) 
{
  int n_thres = tpr.size();
  double p_auc = 0.;
  // Numerical integration of step functions is easy
  for (int j = 1; j < n_thres; j++) {
    p_auc += pauc_tpr_area(tpr, fpr, param, j);
  }
  p_auc = 0.5 * p_auc;
  return p_auc;
}

// double get_perf_pauc_over_fpr(NumericVector &tpr, NumericVector &fpr, NumericVector &param) {
// 
//   int i = 0;
//   double p_auc = 0;
//   // NOTE: first fpr is 1 and param[1] <= 1
//   while (fpr[i] >= param[1]) i++; // fpr starts at 1 and decreases with i
//   // i will be at least 1
//   
// 
//   // check for case that fpr interval considered is very small
//   
//   if (fpr[i] <= param[0]) {
//     p_auc = (tpr[i - 1] - tpr[i]) * (2 - param[1] - param[0]);
//   } 
//   else {
//   
//     // insert right partial area here
//     p_auc = (tpr[i - 1] - tpr[i]) * (2 - param[1] - fpr[i]);
//   
//     while (fpr[i] > param[0]) { // last fpr will be zero and param[0] >= 0 so not out of bounds
//       p_auc += (tpr[i - 1] - tpr[i]) * (2 - fpr[i - 1] - fpr[i]);
//       i++;
//     }
//     
//     p_auc += (tpr[i - 1] - tpr[i]) * (2 - fpr[i - 1] - param[0]);
//   }  
//   
//   p_auc = 0.5 * p_auc;
//   
//   return p_auc;
// }

double get_tpr_at_fixed_fpr(NumericVector &tpr, NumericVector &fpr, NumericVector &param) 
{
  double at = param[0];
  //double out = 0;
  if (at == 1) return param[0];  
  int i = 0;  
  while (fpr[i++] > at);
  //if (fpr[i] == at) out = tpr[i];
  //else
  double out = tpr[i-1];    
  return out;
}

double get_fpr_at_fixed_tpr(NumericVector &tpr, NumericVector &fpr, NumericVector &param) 
{
  double at = param[0];
  if (at == 0) return param[0];
  int i = tpr.size() - 1;  
  while (tpr[i--] < at);
  double out = fpr[i+1];  
  return out;
}
