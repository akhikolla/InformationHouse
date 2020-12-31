#include <Rcpp.h>

using Rcpp::clone;
using Rcpp::NumericVector;
using Rcpp::warning;

// [[Rcpp::export]]
NumericVector QTriC(
    NumericVector p, double min, double max, double mode, bool lower_tail,
    bool log_p)
{
  int n = p.size();

  if (min >= max || mode > max || min > mode)
  {
    warning("NaN(s) produced.");
    return NumericVector(n, R_NaN);
  }

  NumericVector q = clone(p);

  if (log_p)
  {
    q = exp(q);
  }

  if (!lower_tail)
  {
    q = 1.0 - q;
  }

  bool has_nan = false;
  double int_len = max - min;

  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0.0 || q[i] > 1.0)
    {
      q[i] = R_NaN;
      has_nan = true;
    }
    else if (q[i] < (mode - min) / int_len)
    {
      q[i] = min + sqrt(q[i] * int_len * (mode - min));
    }
    else // if (q[i] >= (mode - min) / int_len)
    {
      q[i] = max - sqrt((1.0 - q[i]) * int_len * (max - mode));
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return q;
}

// [[Rcpp::export]]
NumericVector QTriC2(
    NumericVector p, NumericVector min, NumericVector max, NumericVector mode,
    bool lower_tail, bool log_p)
{
  int n = p.size();
  NumericVector q = clone(p);

  if (log_p)
  {
    q = exp(q);
  }

  if (!lower_tail)
  {
    q = 1.0 - q;
  }

  bool has_nan = false;

  for (int i = 0; i < n; i++)
  {
    if (min[i] >= max[i] || mode[i] > max[i] || min[i] > mode[i] ||
        q[i] < 0.0 || q[i] > 1.0)
    {
      q[i] = R_NaN;
      has_nan = true;
    }
    else if (q[i] < (mode[i] - min[i]) / (max[i] - min[i]))
    {
      q[i] = min[i] + sqrt(q[i] * (max[i] - min[i]) * (mode[i] - min[i]));
    }
    else // if (q[i] >= (mode[i] - min[i]) / (max[i] - min[i]))
    {
      q[i] = max[i] - sqrt((1.0 - q[i]) * (max[i] - min[i]) * (max[i] - mode[i]));
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return q;
}
