#' Generating data from a Pareto Distribution.
#'
#' This function is able to generate random Pareto distributed data with 
#' the specified \code{shape} and \code{scale} parameters. The function 
#' has been written to be similar in type to the popular runif and rexp type 
#' of functions for generating data from a particular distribution.
#'
#' @param sample_size number of observations
#' @param shape shape parameter
#' @param scale scale parameter
#'
#' @return Vector of Pareto distributed data of sample size \code{sample_size}
#' with shape parameter \code{shape} and scale parameter \code{scale}.
#' 
#' @examples
#' generate_pareto(10000, 5, 2)
#' generate_pareto(100, 15, 6)
#' 
#' @export
generate_pareto <- function(sample_size, shape, scale){
  U <- runif(sample_size, 0, 1)
  P <- (scale / (U) ^ (1 / shape))
  return (P)
}


#' Obtain estimates for Parameters of Pareto Data from all methods
#'
#' This function combines the results of all the methods (included in this
#' package) provided to estimate the \code{shape} and \code{scale} parameters 
#' of the Pareto data and provides the results in a data frame. Hill's 
#' Estimator is not used in this comparison as it discards a set of 
#' observations. We also note here that when considering the entire data set, 
#' Hill's Estimate is equivalent to the MLE.
#'
#' @param dat vector of observations
#'
#' @return Dataframe with the following columns:
#' \describe{
#'   \item{Method.of.Estimation}{Name of the method used for estimation}
#'   \item{Shape.Parameter}{Estimates of the shape parameter of the data}
#'   \item{Scale.Parameter}{Estimates of the scale parameter of the data}
#' }
#'
#' @examples
#' x <- generate_pareto(10000, 5, 2)
#' generate_all_estimates(x)
#' 
#' @export
generate_all_estimates <- function(dat){

  display <- data.frame("Method of Estimation" = character(),
                        "Shape Parameter" = numeric(),
                        "Scale Parameter" = numeric())

  mle <- alpha_mle(dat)
  ls <- alpha_ls(dat)
  moment <- alpha_moment(dat)
  percentile <- alpha_percentile(dat)
  modified_percentile <- alpha_modified_percentile(dat)
  geometric_percentile <- alpha_geometric_percentile(dat)
  wls <- alpha_wls(dat)

  vector_of_names <- c("Maximum Likelihood Estimate",
                       "Least Squares", "Method of Moments",
                       "Percentiles Method",
                       "Modified Percentiles Method",
                       "Geometric Percentiles Method",
                       "Weighted Least Squares")
  list_of_functions <- list("mle" = mle,
                            "least" = ls, "moment" = moment,
                            "percentile" = percentile,
                            "modified_percentile" = modified_percentile,
                            "geometric_percentiles" = geometric_percentile,
                            "wls" = wls)
  for (i in 1:length(vector_of_names)){
    new_row <- data.frame("Method of Estimation" = vector_of_names[i],
                          "Shape Parameter" = list_of_functions[[i]][[1]],
                          "Scale Parameter" = list_of_functions[[i]][[2]])
    display <- rbind(display, new_row)
  }

  return (display)
}



#' Q-Q Plot to test for Pareto Distribution
#'
#' This function can be used as a first step to identify
#'  whether the data is Pareto distributed before estimating the tail index. If
#'  most of the data points appear to be distributed along a line, it is
#'  possible that the data may be Pareto. Conversely, if most of the data are
#'  distributed non-linearly, then the data is most probably not Pareto
#'  distributed.
#'
#'  This function plots the quantiles of the standard exponential distribution
#'  on the x-axis and the log values of the provided data on the y-axis. If
#'  Pareto data was supplied, a log transformation of this data would result
#'  in an exponential distribution with mean \eqn{\alpha}.
#'  These data points would then show up on the QQ-plot as a line
#'  with slope \eqn{1/\alpha}.
#'
#' The function makes use of the plotly package if available and installed or
#' if not, defaults to the standard R plot.
#'
#' @param dat Data to be tested for Pareto distribution
#'
#' @return A Q-Q plot either using plotly if package is available or else a
#' standard R plot.
#'
#' @examples
#' x <- generate_pareto(10000, 5, 2)
#' pareto_qq_test(x)
#' 
#' @export
pareto_qq_test <- function(dat){

  negative_check(dat)

  sorted_dat <- sort(dat, decreasing = TRUE)
  x_axis <- seq(from = 1, to = length(dat))
  x_axis <- log( (length(dat) + 1) / x_axis)
  y_axis <- log(sorted_dat)

  if ("plotly" %in% installed.packages()[, "Package"] == T){
      #PLOTLY COMMAND
      tmp <- plotly::plot_ly(data.frame(y_axis, x_axis), x = ~ x_axis,
                             y = ~ y_axis, type = "scatter", mode = "markers",
                             marker = list(color = "black"))
      tmp <- plotly::layout(tmp, xaxis = list(title = "theoretical",
                                              zeroline = FALSE),
                     yaxis = list(title = "sample", zeroline = FALSE),
                     title = "<b>QQ Plot</b>")
      plotly::layout(tmp, showlegend = FALSE)
  }
  else{
      #PLOT COMMAND
      plot(x = x_axis, y = y_axis, pch = 19, xlab = "theoretical",
           ylab = "sample", main = "QQ Plot")
      print("Install the plotly package in order to obtain the Q-Q plot with
            more tools.")
  }
}

NULL
