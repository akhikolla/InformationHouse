#' Compute Wasserstein barycenters of subset posteriors
#'
#' This function computes Wasserstein Barycenters of subset posteriors and
#' gives posterior summaries for the full posterior.
#'
#' @param mcmc a three dimensional array (rows =  number of subset posteriors,
#'   columns = number of parameters of the posterior distribution, slices =
#'   samples number of samples for each subset posterior) containing posterior
#'   samples for all subsets
#' @param par.names optional character vector with parameter names
#' @param acc accuracy of the swapping algorithm (default = 0.001)
#' @param iter maximum number of iterations of the swapping algorithm (default = 10)
#' @param out boolean indicating whether output for each iteration of the swapping algorithm should be displayed (default = false)
#'
#' @details The swapping algorithm developed by Puccetti, Rüschendorf and
#'   Vanduffel (2020) is used to compute Wasserstein barycenters of subset
#'   posteriors.
#'
#' @return A \code{wasp} object, which can be further analyzed using the
#'   associated function \code{\link{summary.wasp}}.
#'
#'   A \code{wasp} object contains the following elements (some elements are not
#'   returned if not applicable)
#'
#'   \describe{
#'   \item{\code{barycenter}}{A matrix of posterior samples (rows) for
#'   all parameters (columns) of the full posterior obtained by the swapping algorithm.}
#'   \item{\code{raw}}{An array (\code{dim = c(subsets, parameters, samples)})
#'   containing the raw output from the swapping algorithm.}
#'   \item{\code{call}}{The call to the \code{wasp()} function.}
#'   \item{\code{subsets}}{The amount of subset posteriors in mcmc.}
#'   \item{\code{parameters}}{The amount of parameters in mcmc.}
#'   \item{\code{samples}}{The amount of posterior samples for each subset posterior in mcmc.}
#'   \item{\code{acc}}{Accuracy of the swapping algorithm, default = 0.001.}
#'   \item{\code{iter}}{Maximum amount of iterations for the swapping algorithm, default = 10.}
#'   }
#'
#' @source Puccetti, G., Rüschendorf, L. & Vanduffel, S. (2020). On the
#'   computation of Wasserstein barycenters, Journal of Multivariate Analysis,
#'   176.
#'
#' @examples
#'
#' library(waspr)
#' out <- wasp(pois_logistic,
#'             par.names = c("beta_s", "alpha_l", "beta_l",
#'                           "baseline_sigma", "baseline_mu",
#'                           "correlation", "sigma_s", "sigma_l"))
#' summary(out)
#'
#' @export
#'

wasp <- function(mcmc, par.names = NULL,
                 acc = 0.001, iter = 10, out = FALSE){

  if(!is.null(par.names) & !is.character(par.names)){
    stop("par.names is not a character vector")
  }

  if(!is.numeric(mcmc)){
    stop("mcmc is not numeric")
  }

  if(!is.array(mcmc)){
    stop("mcmc is not a three dimensional array")
  }

  if(length(dim(mcmc)) != 3){
    stop("mcmc is not a three dimensional array")
  }

  subsets = dim(mcmc)[1]
  par = dim(mcmc)[2]
  samples = dim(mcmc)[3]

  if(!is.null(par.names) & par != length(par.names)){
    stop("Names are not provided for each parameter, length(par.names) != dim(mcmc)[2].")
  }

  if(out){
    res_wasp = swap_rcpp(samples = mcmc, acc = acc, iter = iter, out = TRUE)
  }else{
    res_wasp = swap_rcpp(samples = mcmc, acc = acc, iter = iter)
  }


  barycenter = combine(res_wasp)

  if(!is.null(par.names)){
    colnames(barycenter) = par.names
    dimnames(res_wasp) = list(NULL, par.names, NULL)
  }

  call <- match.call()

  output = list('barycenter' = barycenter,
                'raw' = res_wasp,
                'call' = call,
                'subsets' = subsets,
                'parameters' = par,
                'samples' = samples,
                'iter' = iter,
                'acc' = acc)

  class(output) = c("wasp", class(output))

  output
}
