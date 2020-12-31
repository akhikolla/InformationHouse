
#' Get Integerodifferential Quantities
#'
#' Returns an element related to an integral of a derivative of the desired kernel. See references for details.
#'
#' @param X The design matrix
#' @param x The point at which to evaluate.
#' @param The error variance of the model.
#' @param A string, one of "Gaussian", "Matern5_2", and "Matern3_2" indicating the covariance kernel to use.
#' @return The scalar value given the integrated derivative.
#' @references 
#' Mickael Binois, Robert B. Gramacy, Michael Ludkovski (2017), Practical Heteroscedastic Gaussian Process Modeling for Large Simulation Experiments, Journal of Computational and Graphical Statistics
#' @keywords internal
d1 <- function(X, x, sigma, type){
  if(type == "Gaussian"){
    return(d_gauss_cpp(X = X, x = x, sigma = sigma))
  }
  if(type == "Matern5_2"){
    return(d_mat52_cpp(X = X, x = x, sigma = sigma))
  }
  if(type == "Matern3_2"){
    return(d_mat32_cpp(X = X, x = x, sigma = sigma))
  }
}

