#include <Rcpp.h>

using Rcpp::NumericVector;
using Rcpp::warning;

// [[Rcpp::export]]
NumericVector MGTriC(
    NumericVector t, double min, double max, double mode)
{
  int n = t.size();

  if (min >= max || mode >= max || min >= mode)
  {
    warning("NaN(s) produced.");
    return NumericVector(n, R_NaN);
  }

  bool has_nan = false;
  NumericVector m(n);

  for (int i = 0; i < n; i++)
  {
    if (t[i] == 0.0)
    {
      m[i] = R_NaN;
      has_nan = true;
    }
    else
    {
      m[i] = 2.0 * ((max - mode) * exp(min * t[i]) - (max - min) * exp(mode * t[i]) + (mode - min) * exp(max * t[i])) /
             ((max - min) * (mode - min) * (max - mode) *
              pow(t[i], 2));
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return m;
}

// [[Rcpp::export]]
NumericVector MGTriC2(
    NumericVector t, NumericVector min, NumericVector max, NumericVector mode)
{
  int n = t.size();
  bool has_nan = false;
  NumericVector m(n);

  for (int i = 0; i < n; i++)
  {
    if (t[i] == 0.0 ||
        min[i] >= max[i] || mode[i] >= max[i] || min[i] >= mode[i])
    {
      m[i] = R_NaN;
      has_nan = true;
    }
    else
    {
      m[i] = 2.0 * ((max[i] - mode[i]) * exp(min[i] * t[i]) - (max[i] - min[i]) * exp(mode[i] * t[i]) + (mode[i] - min[i]) * exp(max[i] * t[i])) /
             ((max[i] - min[i]) * (mode[i] - min[i]) * (max[i] - mode[i]) *
              pow(t[i], 2));
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return m;
}
