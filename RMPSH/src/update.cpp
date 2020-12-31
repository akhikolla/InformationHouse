#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector update(NumericVector x, double epsilon, double rho, double phi)
{
  int n = x.size();
  int j = 0;
  NumericVector x_updated(2 * n);
  double current_value = 0;
  double updated_value = 0;
  double epsilon_temp = 0;
  float ff = 0;
  int f = 0;

  for (int i = 0; i < (2 * n); ++i)
  {
    j = ceil(float(i / 2));
    epsilon_temp = (1 - 2 * (i % 2)) * epsilon;
    current_value = x[j];
    updated_value = current_value + epsilon_temp;

    if (updated_value > 1 && current_value < 1 - phi)
    {
      ff = log(epsilon_temp / (1 - current_value)) / log(rho);
      f = ceil(ff);
      epsilon_temp = epsilon_temp / pow(rho, f);
      updated_value = current_value + epsilon_temp;
    }
    else if (updated_value < 0 && current_value > phi)
    {
        ff = log(-epsilon_temp / current_value) / log(rho);
        f = ceil(ff);
        epsilon_temp = epsilon_temp / pow(rho, f);
        updated_value = current_value + epsilon_temp;
    }
    if (updated_value > 1 || updated_value < 0)
      updated_value = current_value;
    x_updated[i] = updated_value;
  }

  return x_updated;
}
