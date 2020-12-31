
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% prob %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


###############################################################################@
############################# prob (generic) ###################################
###############################################################################@
#' Calculate the probability of a correct response
#'
#' @description \code{prob} Returns the probability of correct respond to an
#'   item or multiple items with given parameters for a given ability or
#'   abilities, i.e. \eqn{\theta}. For polytomous models, where there are
#'   multiple possible responses, probability of each response category will be
#'   returned.
#'
#' @param ip An \code{\link{Item-class}} or an \code{\link{Itempool-class}}
#'   object containing the item parameters.
#' @param theta An object containing the ability parameters.
#' @param derivative Whether to calculate the first or second derivative of
#'   probability of a response.
#'   \describe{
#'     \item{\code{0}}{No derivative will be calculated. This is the default
#'       value.}
#'     \item{\code{1}}{Calculate the first derivative.}
#'     \item{\code{2}}{Calculate the second derivative.}
#'   }
#' @param expected_value For each possible response value, the probability
#'   of that response is calculated and summed to get the expected value at
#'   a \code{theta} value. Default value is \code{FALSE}.
#'
#' @return Item probabilities at given theta will be returned. If
#'   \code{expected_value} is \code{TRUE}, the expected value(s) of item or
#'   item pool at a given \code{theta} value will be returned.
#'
#' @include Item-class.R
#' @include Itempool-class.R
#'
#' @author Emre Gonulates
#'
setGeneric("prob", function(ip, theta, derivative = 0, expected_value = FALSE)
  {standardGeneric ("prob")})


###############################################################################@
############################# prob (Item) ######################################
###############################################################################@
#' @export
#' @useDynLib irt
#' @rdname prob
#'
#' @examples
#' theta <- rnorm(1)
#' item1 <- generate_item(model = "2PL")
#'
#' # Probability of correct response
#' prob(item1, theta)
#'
#' # First derivative of probability of correct response:
#' prob(item1, theta, derivative = 1)
#'
#' # Second derivative of probability of correct response:
#' prob(item1, theta, derivative = 2)
#'
#' # Probability of each response category for Generalized Partial Credit Model
#' item2 <- generate_item(model = "GPCM", n_categories = 4)
#' prob(item2, theta)
#'
#' # First derivative of each response category
#' prob(item2, theta, derivative = 1)
#'
#' # Second derivative of each response category
#' prob(item2, theta, derivative = 2)
#'
#' # Expected score for a subject with a given theta value
#' prob(item2, theta, expected_value = TRUE)
#'
#' # Probability of each response category for Reparametrized Generalized
#' # Partial Credit Model
#' item3 <- generate_item(model = "GPCM2", n_categories = 3)
#' prob(item3, theta)
#'
#' # Probability of each response category for Graded Response Model
#' item4 <- generate_item(model = "GRM", n_categories = 5)
#' prob(item4, theta)
#'
#' # Multiple theta values
#' theta_n <- rnorm(5)
#'
#' prob(item1, theta_n)
#' prob(item1, theta_n, derivative = 1)
#' prob(item1, theta_n, derivative = 2)
#'
#' prob(item2, theta_n)
#' prob(item2, theta_n, derivative = 1)
#' prob(item2, theta_n, derivative = 2)
#'
setMethod(
  f = "prob", signature = c(ip = "Item"),
  function(ip, theta, derivative = 0, expected_value = FALSE){
    # Expected value can only be calculated when derivative = 0
    if (expected_value && derivative != 0)
      stop("'expected_value' can only be calculated for 'derivative = 0'.")
    # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
    if (ip@model %in%
        names(Pmodels)[sapply(Pmodels, function(x) x$model_family == "UIRT")]) {
      return(prob_4pm_item_cpp(theta = theta, item = ip,
                               derivative = derivative))
    } else if (ip@model %in%
        names(Pmodels)[sapply(Pmodels, function(x) x$model_family == "PIRT")]) {
      result <-  sapply(theta, prob_poly_bare_cpp,  item = ip,
                        derivative = derivative,
                        expected_value = expected_value)
      if (expected_value) {
        result <- c(result)
      } else {
        result <- t(result)
        colnames(result) <- paste0("category.", 0:(ncol(result)-1))
      }
      return(result)
    } else if (ip@model %in% names(Pmodels)[
      sapply(Pmodels, function(x) x$model_family == "MIRT")]) {
      return(prob_mirt_item_cpp(theta = theta, item = ip,
                                derivative = derivative))
    } else stop("This model has not been implemented in 'prob()' function yet.")
  }
)


###############################################################################@
############################# prob (Itempool) #################################
###############################################################################@
#' @export
#'
#' @rdname prob
#'
#' @examples
#'
#' theta <- rnorm(1)
#' ip <- generate_ip(model = "3PL")
#'
#' # Probability of correct response
#' prob(ip, theta)
#'
#' # First derivative of probability of correct response:
#' prob(ip, theta, derivative = 1)
#'
#' # Second derivative of probability of correct response:
#' prob(ip, theta, derivative = 2)
#'
#' # Multiple theta
#' theta_n <- rnorm(5)
#' prob(ip, theta_n)
#' prob(ip, theta_n, derivative = 1)
#' prob(ip, theta_n, derivative = 2)
#'
#'
#' # Probability of each response category for Generalized Partial Credit Model
#' ip <- generate_ip(model = "GPCM", n = 4, n_categories = c(3, 4, 6, 5))
#' prob(ip, theta)
#'
#' # First derivative of each response category
#' prob(ip, theta, derivative = 1)
#'
#' # Second derivative of each response category
#' prob(ip, theta, derivative = 2)
#'
#' # Expected score for a subject with a given theta value for each item
#' prob(ip, theta, expected_value = TRUE)
#'
#' # Probability of a mixture of items models
#' ip <- generate_ip(model = c("GPCM", "2PL", "3PL", "GPCM"),
#'                   n_categories = c(4, 2, 2, 3))
#' prob(ip, theta)
#'
setMethod(
  f = "prob", signature = c(ip = "Itempool"),
  function(ip,  theta, derivative = 0, expected_value = FALSE){
    # Expected value can only be calculated when derivative = 0
    if (expected_value && derivative != 0)
      stop("'expected_value' can only be calculated for 'derivative = 0'.")
    # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
    if (all(ip$model %in%
          names(Pmodels)[sapply(Pmodels,
                                function(x) x$model_family == "UIRT")])) {
      return(prob_4pm_itempool_cpp(theta = theta, ip = ip,
                                   derivative = derivative))
    } else if (all(ip$model %in% names(Pmodels)[
      sapply(Pmodels, function(x) x$model_family == "MIRT")])) {
      return(prob_mirt_itempool_cpp(theta = theta, ip = ip,
                                    derivative = derivative))
    } else if (length(theta) == 1) {
      result <- prob_bare_itempool_cpp(theta, ip, derivative, expected_value)
      dimnames(result) <- list(ip$resp_id, paste0(0:(ncol(result)-1)))
      if (expected_value) result <- result[, 1]
    } else if (expected_value) { # multiple theta and expected_value = TRUE
                                 # mix of dichotomous and polytomous items
      result <- t(sapply(theta, prob_bare_itempool_cpp, ip = ip,
                         derivative = derivative,
                         expected_value = expected_value))
      dimnames(result) <- list(names(theta), ip$resp_id)
    } else
      stop("prob() function cannot be run on item sets that are composed of ",
           "unidimensional and multidimensional items or item sets with ",
           "polytomous items with more than one theta. This model has not ",
           "been implemented in 'prob()' function yet.")
    return(result)
  }
)


###############################################################################@
############################# prob (Testlet) ###################################
###############################################################################@
#' @export
#'
#' @rdname prob
#'
#' @examples
#' theta <- rnorm(1)
#' t1 <- generate_testlet(model_items = "3PL")
#'
#' # Probability of correct response
#' prob(t1, theta)
#'
#' # First derivative of probability of correct response:
#' prob(t1, theta, derivative = 1)
#'
#' # Second derivative of probability of correct response:
#' prob(t1, theta, derivative = 2)
#'
setMethod(
  f = "prob", signature = c(ip = "Testlet"),
  function(ip,  theta, derivative = 0, expected_value = FALSE){
    return(prob(ip = ip@item_list, theta = theta, derivative = derivative,
                expected_value = expected_value))
  }
)


###############################################################################@
############################# prob (REST) ######################################
###############################################################################@
#' @export
#'
#' @rdname prob
setMethod(
  f = "prob", signature = c(ip = "numMatDfListChar"),
  function(ip, theta, derivative = 0, expected_value = FALSE) {
    if (inherits(ip, "numeric")) {
      return(prob(ip = itempool(ip), theta = theta, derivative = derivative,
                  expected_value = expected_value))
    } else if (inherits(ip, c("data.frame", "matrix", "list"))) {
      return(prob(ip = itempool(ip), theta = theta, derivative = derivative,
                  expected_value = expected_value))
    } else {
      stop("Cannot convert object to an 'Item' or an 'Itempool' object. ",
           "Please provide a valid 'Item' or 'Itempool' object using either ",
           "'item()' or 'itempool()' function.")
    }
  }
)

