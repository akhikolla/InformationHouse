/* stable_R.cpp
 *
 * R wrappers
 *
 *
 * Copyright (C) 2015. Javier Royuela del Val
 *                     Federico Simmross Wattenberg
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  Javier Royuela del Val.
 *  E.T.S.I. Telecomunicación
 *  Universidad de Valladolid
 *  Paseo de Belén 15, 47002 Valladolid, Spain.
 *  jroyval@lpi.tel.uva.es
 */

/*
* Check input parameters
*/

#include <Rcpp.h>

using namespace Rcpp;

/* Libstable is writen in C, so we need to make then callable from C++ */
extern "C" {
  #include "stable.h"
}

int checkParams(NumericVector& pars, int parametrization)
{
  int out = 0;

  if (parametrization < 0 || parametrization > 1) {
    perror("Only parametrizations 0 and 1 are accepted");
    out = 5;
  }

  if (pars.size() < 1)
    pars.push_back(2.0); /* alpha */
  else if (pars[0] < 0 || pars[0] > 2.0) {
    perror("Alpha must be between 0.0 and 2.0");
    out = 1;
  }

  if (pars.size() < 2)
    pars.push_back(0.0); /* beta  */
  else if (pars[1] < -1.0 || pars[1] > 1.0) {
    perror("Beta must be between -1.0 and 1.0");
    out = 2;
  }

  if (pars.size() < 3)
    pars.push_back(1.0); /* sigma */
  else if (pars[2] <= 0.0) {
    perror("Sigma must be greater than 0.0");
    out = 3;
  }

  if (pars.size() < 4)
    pars.push_back(0.0); /* mu    */

  return out;
}

/*
* Get stable distribution parameters
*/
NumericVector getPars(StableDist *dist, int parametrization = 0)
{
  NumericVector pars(4);
  pars[0] = dist->alpha;
  pars[1] = dist->beta;
  pars[2] = dist->sigma;
  if (parametrization == 0)
      pars[3] = dist->mu_0;
  else
      pars[3] = dist->mu_1;

  return pars;
}

//' @md
//' @name stable_pdf_and_cdf
//' @aliases stable_cdf stable_pdf
//' @title PDF and CDF of a skew stable distribution.
//' @description Evaluate the PDF or the CDF of the skew stable distribution with parameters
//' pars = c(alpha, beta, sigma, mu) at the points given in x.\cr\cr
//' _parametrization_ argument specifies the parametrization used for the distribution
//' as described by JP Nolan (1997). The default value is _parametrization_ = 0.\cr\cr
//' _tol_ sets the relative error tolerance (precision) to _tol_. The default value is tol = 1e-12.
//' @param x Vector of points where the pdf will be evaluated.
//' @param pars Vector with an initial estimation of the parameters. `pars_init = c(alpha, beta, sigma, mu)`, where
//' * alpha: shape / stability parameter, with 0 < alpha <= 2.
//' * beta: skewness parameter, with -1 <= beta <= 1.
//' * sigma: scale parameter, with 0 < sigma.
//' * mu: location parameter, with mu real.
//' @param parametrization Parametrization used for the skew stable distribution, as defined by JP Nolan (1997). By default, parametrization = 0.
//' @param tol Relative error tolerance (precission) of the calculated values. By default, tol = 1e-12.
//' @return A numeric vector.
//' @author Javier Royuela del Val, Federico Simmross Wattenberg and Carlos Alberola López\cr\cr
//'         Maintainer: Javier Royuela del Val <jroyval@@lpi.tel.uva.es>
//' @references Nolan JP (1997). Numerical Calculation of Stable Densities and Distribution Functions. Stochastic Models, 13(4) 759-774.
//' @keywords distribution
//' @examples
//' pars <- c(1.5, 0.9, 1, 0)
//' x <- seq(-5, 10, 0.001)
//'
//' pdf <- stable_pdf(x, pars)
//' cdf <- stable_cdf(x, pars)
//'
//' plot(x, pdf, type = "l")
//' @export
// [[Rcpp::export]]
NumericVector stable_pdf(NumericVector x, NumericVector pars, int parametrization=0, double tol=1e-12)
{
  NumericVector out(x.size());

  if(checkParams(pars, parametrization) != 0) {
    perror("No valid parameters provided");
    out.fill(NA_REAL);
    return out;
  }

  StableDist * dist = stable_create(pars[0], pars[1], pars[2], pars[3], parametrization);
  stable_set_relTOL(tol);
  stable_pdf(dist, &(x[0]), x.size(), &(out[0]), NULL);
  stable_free(dist);

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector stable_cdf(NumericVector x, NumericVector pars, int parametrization=0, double tol=1e-12)
{
  NumericVector out(x.size());

  if(checkParams(pars, parametrization) != 0) {
    perror("No valid parameters provided");
    out.fill(NA_REAL);
    return out;
  }

  StableDist * dist = stable_create(pars[0], pars[1], pars[2], pars[3], parametrization);
  stable_set_relTOL(tol);
  stable_cdf(dist, &(x[0]), x.size(), &(out[0]), NULL);
  stable_free(dist);

  return out;
}

//' @md
//' @title Quantile function of skew stable distributions
//' @description Evaluate the quantile function (CDF^-1) of the skew stable distribution
//' with parameters pars = c(alpha, beta, sigma, mu) at the points given in p.\cr\cr
//' _parametrization_ argument specifies the parametrization used for the distribution
//' as described by JP Nolan (1997). The default value is _parametrization_ = 0.\cr\cr
//' _tol_ sets the relative error tolerance (precission) to _tol_. The default value is tol = 1e-12.
//' @param p Vector of points where the quantile function will be evaluated, with  0 < p\[i] < 1.0
//' @param pars Vector with an initial estimation of the parameters. `pars_init = c(alpha, beta, sigma, mu)`, where
//' * alpha: shape / stability parameter, with 0 < alpha <= 2.
//' * beta: skewness parameter, with -1 <= beta <= 1.
//' * sigma: scale parameter, with 0 < sigma.
//' * mu: location parameter, with mu real.
//' @param parametrization Parametrization used for the skew stable distribution, as defined by JP Nolan (1997). By default, parametrization = 0.
//' @param tol Relative error tolerance (precission) of the calculated values. By default, tol = 1e-12.
//' @return A numeric vector.
//' @author Javier Royuela del Val, Federico Simmross Wattenberg and Carlos Alberola López\cr\cr
//'         Maintainer: Javier Royuela del Val <jroyval@@lpi.tel.uva.es>
//' @keywords distribution
//' @export
// [[Rcpp::export]]
NumericVector stable_q(NumericVector p, NumericVector pars, int parametrization=0, double tol=1e-12)
{
  NumericVector out(p.size());

  if(checkParams(pars, parametrization) != 0) {
    perror("No valid parameters provided");
    out.fill(NA_REAL);
    return out;
  }

  int k = 0;
  int N = p.size();
  for (k = 0; k < N; k++) {
    if (p[k]>1 || p[k]<0) {
      perror("p but must be between 0 and 1");
      out.fill(NA_REAL);
      return out;
    }
  }

  StableDist * dist = stable_create(pars[0], pars[1], pars[2], pars[3], parametrization);
  stable_set_relTOL(tol);
  stable_q(dist, &(p[0]), p.size(), &(out[0]), NULL);
  stable_free(dist);

  return out;
}

//' @md
//' @title Skew stable distribution random sample generation.
//' @description `stable_rnd(N, pars)` generates N random samples of a skew stable distribuiton
//' with parameters pars = c(alpha, beta, sigma, mu) using the Chambers, Mallows,
//' and Stuck (1976) method.\cr\cr
//' @param N Number of values to generate.
//' @param pars Vector with an initial estimation of the parameters. `pars_init = c(alpha, beta, sigma, mu)`, where
//' * alpha: shape / stability parameter, with 0 < alpha <= 2.
//' * beta: skewness parameter, with -1 <= beta <= 1.
//' * sigma: scale parameter, with 0 < sigma.
//' * mu: location parameter, with mu real.
//' @param parametrization Parametrization used for the skew stable distribution, as defined by JP Nolan (1997). By default, parametrization = 0.
//' @references Chambers JM, Mallows CL, Stuck BW (1976). A Method for Simulating Stable Random Variables. Journal of the American Statistical Association, 71(354), 340-344. doi:10.1080/01621459.1976.10480344.
//' @return A numeric vector.
//' @author Javier Royuela del Val, Federico Simmross Wattenberg and Carlos Alberola López\cr\cr
//'         Maintainer: Javier Royuela del Val <jroyval@@lpi.tel.uva.es>
//' @keywords distribution
//' @examples
//' N <- 1000
//' pars <- c(1.25, 0.95, 1.0, 0.0)
//' rnd <- stable_rnd(N, pars)
//'
//' hist(rnd)
//' @export
// [[Rcpp::export]]
NumericVector stable_rnd(int N, NumericVector pars, int parametrization=0)
{
  NumericVector out(N);

  if(checkParams(pars, parametrization) != 0) {
    perror("No valid parameters provided");
    out.fill(NA_REAL);
    return out;
  }

  StableDist * dist = stable_create(pars[0], pars[1], pars[2], pars[3], parametrization);
  stable_rnd(dist, &(out[0]), N);
  stable_free(dist);

  return out;
}

//' @md
//' @title Methods for parameter estimation of skew stable distributions.
//' @name stable_fit
//' @aliases stable_fit_init stable_fit_koutrouvelis stable_fit_mle stable_fit_mle2d
//' @description A set of functions are provided that perform the parameter estimation of skew stable distributions with different methods.
//' @details
//' * `stable_fit_init()` uses McCulloch's method of quantiles \[3]. This is usually a good initialization for the rest of the methods.
//' * `stable_fit_koutrouvelis()` implements Koutrouvellis' method based on the characteristic function \[4].
//' * `stable_fit_mle()` implements a Maximum likelihood estimation.
//' * `stable_fit_mle2()` implements a modified maximum likelihood estimation as described in \[1].
//' @param rnd Random sample
//' @param parametrization Parametrization used for the skew stable distribution, as defined by JP Nolan (1997). By default, parametrization = 0.
//' @return A numeric vector.
//' @references
//' * \[1] Royuela-del-Val J, Simmross-Wattenberg F, Alberola López C (2017). libstable: Fast, Parallel and High-Precision Computation of alpha-stable Distributions in R, C/C++ and MATLAB. Journal of Statistical Software, 78(1), 1-25. doi:10.18637/jss.v078.i01
//' * \[2] Chambers JM, Mallows CL, Stuck BW (1976). A Method for Simulating Stable Random Variables. Journal of the American Statistical Association, 71(354), 340-344. doi:10.1080/01621459.1976.10480344.
//' * \[3] McCulloch JH (1986). Simple Consistent Estimators of Stable Distribution Parameters. Communications in Statistics - Simulation and Computation, 15(4), 1109-1136. doi:10.1080/03610918608812563.
//' * \[4] Koutrouvelis IA (1981). An Iterative Procedure for the Estimation of the Parameters of Stable Laws. Communications in Statistics - Simulation and Computation, 10(1), 17-28. doi:10.1080/03610918108812189.
//' * \[5] Nolan JP (1997). Numerical Calculation of Stable Densities and Distribution Functions. Stochastic Models, 13(4) 759-774. doi:10.1080/15326349708807450.
//' @keywords distribution
//' @author Javier Royuela del Val, Federico Simmross Wattenberg and Carlos Alberola López\cr\cr
//'         Maintainer: Javier Royuela del Val <jroyval@@lpi.tel.uva.es>
//' @examples
//' # Set alpha, beta, sigma and mu stable parameters in a vector
//' pars <- c(1.5, 0.9, 1, 0)
//'
//' # Generate 300 random values
//' rnd <- stable_rnd(300, pars)
//'
//' # Estimate the parameters of the skew stable distribution given
//' # the generated sample:
//'
//' # Using the McCulloch's estimator:
//' pars_init <- stable_fit_init(rnd)
//'
//' # Using the Koutrouvelis' estimator, with McCulloch estimation
//' # as a starting point:
//' pars_est_K <- stable_fit_koutrouvelis(rnd, pars_init)
//'
//' # Using maximum likelihood estimator:
//' # pars_est_ML <- stable_fit_mle(rnd, pars_est_K)
//'
//' # Using modified maximum likelihood estimator (see [1]):
//' # pars_est_ML2 <- stable_fit_mle2d(rnd, pars_est_K)
//' @export
// [[Rcpp::export]]
NumericVector stable_fit_init(NumericVector rnd, int parametrization=0)
{
  /* Non an iterative method. Parameters do not influence the result */
  StableDist* dist = stable_create(2.0, 0.0, 1.0, 0.0, 0);
  stable_fit_init(dist, &(rnd[0]), rnd.size(), NULL, NULL);

  NumericVector out = getPars(dist, parametrization);

  stable_free(dist);
  return out;
}

//' @md
//' @rdname stable_fit
//' @param pars_init Vector with an initial estimation of the parameters. `pars_init = c(alpha, beta, sigma, mu)`, where
//' * alpha: shape / stability parameter, with 0 < alpha <= 2.
//' * beta: skewness parameter, with -1 <= beta <= 1.
//' * sigma: scale parameter, with 0 < sigma.
//' * mu: location parameter, with mu real.
//' @export
// [[Rcpp::export]]
NumericVector stable_fit_koutrouvelis(NumericVector rnd, NumericVector pars_init = NumericVector::create(), int parametrization=0)
{
  if (pars_init.size() == 0)
    pars_init = stable_fit_init(rnd, parametrization);

  if(checkParams(pars_init, parametrization) != 0) {
    perror("No valid parameters provided");
    NumericVector out(4, NA_REAL);
    return out;
  }

  StableDist* dist = stable_create(pars_init[0], pars_init[1], pars_init[2], pars_init[3], parametrization);


  if (stable_fit_koutrouvelis(dist, &(rnd[0]), rnd.size()) < 0) {
      Rprintf("Stable_fit_koutrouvelis error");
  }

  NumericVector out = getPars(dist, parametrization);

  stable_free(dist);
  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector stable_fit_mle(NumericVector rnd, NumericVector pars_init = NumericVector::create(), int parametrization = 0)
{

  if (pars_init.size() == 0) {
    pars_init = stable_fit_init(rnd, parametrization);
    Rprintf("INIT MCCULLCOH\n");
  }
  else {
    Rprintf("SKIP INIT\n");
  }

  if(checkParams(pars_init, parametrization) != 0) {
    perror("No valid parameters provided");
    NumericVector out(4, NA_REAL);
    return out;
  }

  StableDist* dist = stable_create(pars_init[0], pars_init[1], pars_init[2], pars_init[3], parametrization);

  if (stable_fit_mle(dist, &(rnd[0]), rnd.size()) < 0) {
      Rprintf("Stable_fit_mle error");
  }

  NumericVector out = getPars(dist, parametrization);

  stable_free(dist);
  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector stable_fit_mle2d(NumericVector rnd, NumericVector pars_init = NumericVector::create(), int parametrization = 0)
{

  if (pars_init.size() == 0)
    pars_init = stable_fit_init(rnd, parametrization);

  if(checkParams(pars_init, parametrization) != 0) {
    perror("No valid parameters provided");
    NumericVector out(4, NA_REAL);
    return out;
  }

  StableDist* dist = stable_create(pars_init[0], pars_init[1], pars_init[2], pars_init[3], parametrization);
  stable_fit_mle2d(dist, &(rnd[0]), rnd.size());

  NumericVector out = getPars(dist, parametrization);

  stable_free(dist);
  return out;
}
