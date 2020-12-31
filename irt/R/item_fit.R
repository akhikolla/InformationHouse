

###############################################################################@
############################# item_fit #########################################
###############################################################################@
#' Calculate item-fit indices
#' @description
#' \code{item_fit} calculates the fit of an item to a given psychometric model.
#'
#'
#' @param ip An \code{\link{Itempool-class}} object.
#' @param resp A vector of item responses.
#' @param theta An vector containing ability parameters.
#' @param type The type of the item-fit index.
#'
#' @return A vector of item-fit index values.
#'
#' @include Item-class.R
#' @include Itempool-class.R
#' @include Item-class-methods.R
#' @include Itempool-class-methods.R
#'
#' @author Emre Gonulates
#'
#' @export
#'
item_fit <- function(ip, resp, theta, type = "Q3") {
  if (type == "Q3") {
    d <- resp - prob(ip = ip, theta = theta)
    return(stats::cor(d, use = "pairwise.complete.obs"))
  }
}
