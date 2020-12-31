#include <Rcpp.h>
#include <dqrng.h>

using dqrng::dqrunif;
using Rcpp::NumericVector;
using Rcpp::runif;
using Rcpp::warning;

// [[Rcpp::export]]
NumericVector RTriC(int n, double min, double max, double mode, bool dqrng)
{
  if (min >= max || mode > max || min > mode)
  {
    warning("NaN(s) produced.");
    return NumericVector(n, R_NaN);
  }

  NumericVector r = dqrng ? dqrunif(n) : runif(n);
  double int_len = max - min;

  for (int i = 0; i < n; i++)
  {
    if (r[i] < (mode - min) / int_len)
    {
      r[i] = min + sqrt(r[i] * int_len * (mode - min));
    }
    else // if (r[i] >= (mode - min) / int_len)
    {
      r[i] = max - sqrt((1.0 - r[i]) * int_len * (max - mode));
    }
  }

  return r;
}

// [[Rcpp::export]]
NumericVector RTriC2(
    int n, NumericVector min, NumericVector max, NumericVector mode, bool dqrng)
{
  NumericVector r = dqrng ? dqrunif(n) : runif(n);
  bool has_nan = false;

  for (int i = 0; i < n; i++)
  {
    if (min[i] >= max[i] || mode[i] > max[i] || min[i] > mode[i])
    {
      r[i] = R_NaN;
      has_nan = true;
    }
    else if (r[i] < (mode[i] - min[i]) / (max[i] - min[i]))
    {
      r[i] = min[i] + sqrt(r[i] * (max[i] - min[i]) * (mode[i] - min[i]));
    }
    else // if (r[i] >= (mode[i] - min[i]) / (max[i] - min[i]))
    {
      r[i] = max[i] - sqrt((1.0 - r[i]) * (max[i] - min[i]) *
                           (max[i] - mode[i]));
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return r;
}
