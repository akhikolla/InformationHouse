#include <Rcpp.h>
using namespace Rcpp;

//' Interpolation-specific sequencer
//'
//' @param samples NumericVector. A data stream.
//' @keywords internal
// [[Rcpp::export]]
NumericVector zero2one(NumericVector samples) {
  double length_out = samples.size();
  double by = 1/(length_out - 1);
  NumericVector result = 0;
  for (int i= 0; i < length_out; ++i) {
    double new_result = i * by;
    result.push_back(new_result);
  }
  return result;
}

//' Match a resampled interval proportion to its corresponding
//' originally-sampled interval in C++.
//'
//' @param proportion double. The interval proportion to match.
//' @param reference_frame DataFrame containing a \code{prop_min} column to
//'   serve as originally-sampled interval proportions
//'
//' @keywords internal
// [[Rcpp::export]]
int interval_match(
    double proportion, NumericVector references
) {
  int index = references.size() - 1;
  if (proportion > references[index]) return index;
  while ((proportion <= references[index]) & (index >= 0)) {
    index--;
  }
  if (index < 0) ++index;
  return index;
}

//' @rdname sensor_resample
//' @aliases Resample Interpolate
//' @keywords internal
// [[Rcpp::export]]
NumericVector interpolate_C(
    NumericVector original_samples, int target_frequency
) {

  if (original_samples.size() == target_frequency) {
    return original_samples;
  }

  NumericVector intervals = zero2one(original_samples);

  // Interpolation information

    NumericVector prop_min(0);
    for (int i=0; i < (intervals.size() - 1); ++i) {
      prop_min.push_back(intervals[i]);
    }

    NumericVector start(0);
    for (int i=0; i < (original_samples.size() - 1); ++i) {
      start.push_back(original_samples[i]);
    }

    NumericVector rise = diff(original_samples);

    NumericVector run = diff(intervals);

  // New Data

    NumericVector proportion(0);
    for (double i = 0; i < target_frequency; ++i) {
      double new_result = i / target_frequency;
      proportion.push_back(new_result);
    }

    IntegerVector index(0);
    for (int i = 0; i < proportion.size(); ++i) {
      int new_index = interval_match(
        proportion[i], prop_min
      );
      index.push_back(new_index);
    }

    NumericVector window_fraction(0);
    for (int i = 0; i < index.size(); ++i) {
      int j = index[i];
      double new_fraction = (proportion[i] - prop_min[j])/run[j];
      window_fraction.push_back(new_fraction);
    }

    if(is_true(any(window_fraction < 0))) {
      stop("Detected window fractions < 0");
    }
    if(is_true(any(window_fraction > 1))) {
      stop("Detected window fractions > 1");
    }

    NumericVector final_values(0);
    for (int i = 0; i < window_fraction.size(); ++i) {
      int j = index[i];
      double new_value =  start[j] + (
        rise[j] * window_fraction[i]
      );
      final_values.push_back(new_value);
    }

  return final_values;

}
