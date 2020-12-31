#' LibstableR: Fast and accurate evaluation, random number generation
#' and parameter estimation of skew stable distributions.
#'
#' LibstableR provides functions to work with skew stable distributions
#' in a fast and accurate way \[1]. It performs:
#'
#' * Fast and accurate evaluation of the probability density function (PDF) and cumulative density function (CDF).
#' * Fast and accurate evaluation of the quantile function (CDF^{-1}).
#' * Random numbers generation \[2].
#' * Skew stable parameter estimation with:
#'     * McCulloch's method of quantiles \[3].
#'     * Koutrouvellis' method based on the characteristic function \[4].
#'     * Maximum likelihood estimation.
#'     * Modified maximum likelihood estimation as described in \[1].
#' *The evaluation of the PDF and CDF is based on the formulas provided by John P Nolan in \[5].
#'
#' @md
#' @author Javier Royuela del Val, Federico Simmross Wattenberg and Carlos Alberola López;\cr\cr
#'         Maintainer: Javier Royuela del Val <jroyval@@lpi.tel.uva.es>
#' @references
#' * \[1] Royuela-del-Val J, Simmross-Wattenberg F, Alberola López C (2017). libstable: Fast, Parallel and High-Precision Computation of alpha-stable Distributions in R, C/C++ and MATLAB. Journal of Statistical Software, 78(1), 1-25. doi:10.18637/jss.v078.i01
#' * \[2] Chambers JM, Mallows CL, Stuck BW (1976). A Method for Simulating Stable Random Variables. Journal of the American Statistical Association, 71(354), 340-344. doi:10.1080/01621459.1976.10480344
#' * \[3] McCulloch JH (1986). Simple Consistent Estimators of Stable Distribution Parameters. Communications in Statistics - Simulation and Computation, 15(4), 1109-1136. doi:10.1080/03610918608812563
#' * \[4] Koutrouvelis IA (1981). An Iterative Procedure for the Estimation of the Parameters of Stable Laws. Communications in Statistics - Simulation and Computation, 10(1), 17-28. doi:10.1080/03610918108812189
#' * \[5] Nolan JP (1997). Numerical Calculation of Stable Densities and Distribution Functions. Stochastic Models, 13(4), 759-774. doi:10.1080/15326349708807450
#' @name libstableR
#' @docType package
#' @keywords package
#' @useDynLib libstableR, .registration=TRUE
#' @importFrom Rcpp sourceCpp evalCpp
#' @examples
#' # Set alpha, beta, sigma and mu stable parameters in a vector
#' pars <- c(1.5, 0.9, 1, 0)
#'
#' # Generate an abscissas axis and probabilities vector
#' x <- seq(-5, 10, 0.05)
#' p <- seq(0.01, 0.99, 0.01)
#'
#' # Calculate pdf, cdf and quantiles
#' pdf <- stable_pdf(x, pars)
#' cdf <- stable_cdf(x, pars)
#' xq  <- stable_q(p, pars)
#'
#' # Generate 300 random values
#' rnd <- stable_rnd(300, pars)
#'
#' # Estimate the parameters of the skew stable distribution given
#' # the generated sample:
#'
#' # Using the McCulloch's estimator:
#' pars_est_M <- stable_fit_init(rnd)
#'
#' # Using the Koutrouvelis' estimator:
#' pars_est_K <- stable_fit_koutrouvelis(rnd, pars_est_M)
#'
#' # Using maximum likelihood estimator, with McCulloch estimation
#' # as a starting point:
#' # pars_est_ML <- stable_fit_mle(rnd, pars_est_M)
#'
#' # Using modified maximum likelihood estimator (See [1]):
#' # pars_est_ML2 <- stable_fit_mle2d(rnd, pars_est_M)
NULL
