#include <Rcpp.h>

using Rcpp::ComplexVector;
using Rcpp::NumericVector;
using Rcpp::warning;
using std::complex;

// [[Rcpp::export]]
ComplexVector CTriC(NumericVector t, double min, double max, double mode)
{
  int n = t.size();

  if (min >= max || mode >= max || min >= mode)
  {
    warning("NaN(s) produced.");
    return ComplexVector(n, complex<double>(R_NaN, 0.0));
  }

  bool has_nan = false;
  complex<double> x(0.0, 1.0);
  ComplexVector c(n);
  Rcomplex rc;

  for (int i = 0; i < n; i++)
  {
    if (t[i] == 0.0)
    {
      rc.r = R_NaN;
      rc.i = 0.0;
      c[i] = rc;
      has_nan = true;
    }
    else
    {
      complex<double> cc = -2.0 * ((max - mode) * exp(x * min * t[i]) - (max - min) * exp(x * mode * t[i]) + (mode - min) * exp(x * max * t[i])) /
                           ((max - min) * (mode - min) * (max - mode) * pow(t[i], 2));
      rc.r = cc.real();
      rc.i = cc.imag();
      c[i] = rc;
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return c;
}

// [[Rcpp::export]]
ComplexVector CTriC2(
    NumericVector t, NumericVector min, NumericVector max, NumericVector mode)
{
  int n = t.size();
  bool has_nan = false;
  complex<double> x(0.0, 1.0);
  ComplexVector c(n);
  Rcomplex rc;

  for (int i = 0; i < n; i++)
  {
    if (t[i] == 0.0 || min[i] >= max[i] || mode[i] >= max[i] ||
        min[i] >= mode[i])
    {
      rc.r = R_NaN;
      rc.i = 0.0;
      c[i] = rc;
      has_nan = true;
    }
    else
    {
      complex<double> cc = -2.0 * ((max[i] - mode[i]) * exp(x * min[i] * t[i]) - (max[i] - min[i]) * exp(x * mode[i] * t[i]) + (mode[i] - min[i]) * exp(x * max[i] * t[i])) /
                           ((max[i] - min[i]) * (mode[i] - min[i]) * (max[i] - mode[i]) *
                            pow(t[i], 2));
      rc.r = cc.real();
      rc.i = cc.imag();
      c[i] = rc;
    }
  }

  if (has_nan)
  {
    warning("NaN(s) produced.");
  }

  return c;
}
