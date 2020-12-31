#include <Rcpp.h>

using Rcpp::NumericVector;
using Rcpp::warning;

// [[Rcpp::export]]
NumericVector ESTriC(
    NumericVector p, double min, double max, double mode, bool lower_tail,
    bool log_p)
{
  int n = p.size();

  if (min >= max || mode > max || min > mode)
  {
    warning("NaN(s) produced.");
    return NumericVector(n, R_NaN);
  }

  bool has_nan = false;
  NumericVector es(n);

  for (int i = 0; i < n; i++)
  {
    if (log_p)
    {
      p[i] = exp(p[i]);
    }

    if (!lower_tail)
    {
      p[i] = 1.0 - p[i];
    }

    if (p[i] == 0.0 || p[i] < 0.0 || p[i] > 1.0)
    {
      es[i] = R_NaN;
      has_nan = true;
    }
    else if (p[i] < (mode - min) / (max - min))
    {
      es[i] = ((p[i] * min) + (2.0 / 3.0) * sqrt((max - min) * (mode - min)) *
                                  pow(p[i], 1.5)) /
              p[i];
    }
    else // if (p[i] >= (mode - min) / (max - min))
    {
      double b = (mode - min) / (max - min);
      es[i] = ((b * min) + (2.0 / 3.0) * sqrt((max - min) * (mode - min)) * pow(b, 1.5) + (((p[i] * max) + (2.0 / 3.0) * sqrt((max - min) * (max - mode)) * pow((1.0 - p[i]), 1.5)) - ((b * max) + (2.0 / 3.0) * sqrt((max - min) * (max - mode)) * pow((1.0 - b), 1.5)))) / p[i];
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return es;
}

// [[Rcpp::export]]
NumericVector ESTriC2(
    NumericVector p, NumericVector min, NumericVector max, NumericVector mode,
    bool lower_tail, bool log_p)
{
  int n = p.size();
  bool has_nan = false;
  NumericVector es(n);

  for (int i = 0; i < n; i++)
  {
    if (log_p)
    {
      p[i] = exp(p[i]);
    }

    if (!lower_tail)
    {
      p[i] = 1.0 - p[i];
    }

    if (min[i] >= max[i] || mode[i] > max[i] || min[i] > mode[i] ||
        p[i] <= 0.0 || p[i] > 1.0)
    {
      es[i] = R_NaN;
      has_nan = true;
    }
    else if (p[i] < (mode[i] - min[i]) / (max[i] - min[i]))
    {
      es[i] = ((p[i] * min[i]) + (2.0 / 3.0) * sqrt((max[i] - min[i]) * (mode[i] - min[i])) * pow(p[i], 1.5)) / p[i];
    }
    else // if (p[i] >= (mode - min) / (max - min))
    {
      double b = (mode[i] - min[i]) / (max[i] - min[i]);
      es[i] = ((b * min[i]) + (2.0 / 3.0) * sqrt((max[i] - min[i]) * (mode[i] - min[i])) * pow(b, 1.5) + (((p[i] * max[i]) + (2.0 / 3.0) * sqrt((max[i] - min[i]) * (max[i] - mode[i])) * pow((1.0 - p[i]), 1.5)) - ((b * max[i]) + (2.0 / 3.0) * sqrt((max[i] - min[i]) * (max[i] - mode[i])) * pow((1.0 - b), 1.5)))) / p[i];
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return es;
}
