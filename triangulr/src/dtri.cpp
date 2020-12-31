#include <Rcpp.h>

using Rcpp::NumericVector;
using Rcpp::warning;

// [[Rcpp::export]]
NumericVector DTriC(
    NumericVector x, double min, double max, double mode, bool log)
{
  int n = x.size();

  if (min >= max || mode > max || min > mode)
  {
    warning("NaN(s) produced.");
    return NumericVector(n, R_NaN);
  }

  NumericVector d(n);

  for (int i = 0; i < n; i++)
  {
    if (x[i] < min || x[i] > max)
    {
      d[i] = 0.0;
    }
    else if (min <= x[i] && x[i] < mode)
    {
      d[i] = 2.0 * (x[i] - min) / ((max - min) * (mode - min));
    }
    else if (x[i] == mode)
    {
      d[i] = 2.0 / (max - min);
    }
    else // if (mode < x[i] && x[i] <= max)
    {
      d[i] = 2.0 * (max - x[i]) / ((max - min) * (max - mode));
    }
  }

  if (log)
  {
    return Rcpp::log(d);
  }

  return d;
}

// [[Rcpp::export]]
NumericVector DTriC2(
    NumericVector x, NumericVector min, NumericVector max, NumericVector mode,
    bool log)
{
  int n = x.size();
  bool has_nan = false;
  NumericVector d(n);

  for (int i = 0; i < n; i++)
  {
    if (min[i] >= max[i] || mode[i] > max[i] || min[i] > mode[i])
    {
      d[i] = R_NaN;
      has_nan = true;
    }
    else if (x[i] < min[i] || x[i] > max[i])
    {
      d[i] = 0.0;
    }
    else if (min[i] <= x[i] && x[i] < mode[i])
    {
      d[i] = 2.0 * (x[i] - min[i]) / ((max[i] - min[i]) * (mode[i] - min[i]));
    }
    else if (x[i] == mode[i])
    {
      d[i] = 2.0 / (max[i] - min[i]);
    }
    else // if (mode[i] < x[i] && x[i] <= max[i])
    {
      d[i] = 2.0 * (max[i] - x[i]) / ((max[i] - min[i]) * (max[i] - mode[i]));
    }
  }

  if (log)
  {
    return Rcpp::log(d);
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return d;
}
