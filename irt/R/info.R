

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


###############################################################################@
############################# info (generic) ###################################
###############################################################################@
#' Calculates the information of an "Item" object
#' @description
#' This function sets a generic method for calculating the information of
#' a suitable object
#'
#' @param ip An \code{\link{Item-class}}, \code{\link{Itempool-class}} or
#'   \code{\link{Testlet-class}} object.
#' @param theta An vector of ability parameters.
#' @param tif If it is \code{TRUE}, function will return total
#' information obtained from each item for a given theta. It simply adds
#' information of individual items.
#' @param observed If \code{TRUE}, observed information calculated instead of
#'   the default expected information.
#' @param resp A response string (vector or a matrix). Necessary for observed
#'   information.
#'
#' @return A vector (or matrix) consist of item or test information.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#' @include Itempool-class-methods.R
#'
#' @author Emre Gonulates
#'
#'
setGeneric("info", function(ip, theta, tif = FALSE, observed = FALSE,
                            resp = NULL)
  {standardGeneric ("info")})


###############################################################################@
############################# info (Item) ######################################
###############################################################################@
#' @export
#'
#' @rdname info
#'
#' @examples
#' info(ip = generate_item(model = "Rasch"), theta = rnorm(1))
#' info(ip = generate_item(model = "1PL"), theta = rnorm(1))
#' info(ip = generate_item(model = "2PL"), theta = rnorm(1))
#' info(ip = generate_item(model = "3PL"), theta = rnorm(1))
#' info(ip = generate_item(model = "4PL"), theta = rnorm(1))
#' info(ip = generate_item(model = "GRM"), theta = rnorm(1))
#' info(ip = generate_item(model = "GPCM"), theta = rnorm(1))
#' info(ip = generate_item(model = "PCM"), theta = rnorm(1))
#' info(ip = generate_item(model = "GPCM2"), theta = rnorm(1))
#'
setMethod(
  f = "info", signature = c(ip = "Item"),
  function(ip, theta, tif = FALSE, observed = FALSE, resp = NULL){
    # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
    if (all(ip$model %in%
              c(names(Pmodels)[
                sapply(Pmodels, function(x) x$model_family == 'UIRT')],
                "GRM", "GPCM", "PCM", "GPCM2"))) {
      return(info_item_cpp(theta = theta, item = ip, observed = observed,
                           resp = resp))
    } else stop("This model has not been implemented in 'info()' function.",
                call. = FALSE)
  }
)


############################# .info_itempool ##############################@###
.info_itempool <-  function(ip, theta, tif = FALSE, observed = FALSE,
                             resp = NULL) {
  if (inherits(ip, "Testlet")) ip <- ip@item_list
  # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
  acceptable_models <- c(
    names(Pmodels)[sapply(Pmodels, function(x) x$model_family == 'UIRT')],
    "GRM", "GPCM", "PCM", "GPCM2")
  # All items should be acceptable
  test <- FALSE

  # All standalone Item models should be valid
  if ("Item" %in% sapply(ip, class))
    test <- test | all(ip[sapply(ip, class)  == "Item"]$model %in%
                         acceptable_models)
  # Check whether Testlet items are valid
  if ("Testlet" %in% sapply(ip, class))
    test <- test | all(
      unlist(lapply(ip[sapply(ip, class)  == "Testlet"],
                    function(x) x@item_list$model)) %in% acceptable_models)
  # Check whether resp is acceptable
  if (!is.null(resp)) {
    if (!is.numeric(resp))
      stop("Invalid response vector. It should be numeric.")
    if (!is.matrix(resp)) resp <- matrix(resp, nrow = 1)
    ip_size <- get_itempool_size(ip)
    if (nrow(resp) != length(theta) || ncol(resp) != ip_size["items"])
      stop("Invalid response vector/matrix. The item pool size and response ",
           "size do not match.")
  }
  if (test) {
    return(info_itempool_cpp(theta = theta, ip = ip, tif = tif,
                              observed = observed, resp = resp))
  } else stop("This model is not implemented in 'info()' function.")
}


###############################################################################@
############################# info (Itempool) #################################
###############################################################################@
#' @export
#' @rdname info
#'
#' @examples
#' info(ip = generate_ip(model = "Rasch"), theta = rnorm(1))
#' info(ip = generate_ip(model = "1PL"), theta = rnorm(1))
#' info(ip = generate_ip(model = "2PL"), theta = rnorm(1))
#' info(ip = generate_ip(model = "3PL"), theta = rnorm(1))
#' info(ip = generate_ip(model = "4PL"), theta = rnorm(1))
#' info(ip = generate_ip(model = "GRM"), theta = rnorm(1))
#' info(ip = generate_ip(model = "GPCM"), theta = rnorm(1))
#' info(ip = generate_ip(model = "PCM"), theta = rnorm(1))
#' info(ip = generate_ip(model = "GPCM2"), theta = rnorm(1))
#'
#' # Multiple Thetas
#' info(ip = generate_ip(model = "3PL"), theta = rnorm(5))
#' info(ip = generate_ip(model = "GRM"), theta = rnorm(7))
#'
#' # Test information function value at theta
#' info(ip = generate_ip(model = "3PL"), theta = rnorm(5), tif = TRUE)
#' info(ip = generate_ip(model = "GRM"), theta = rnorm(7), tif = TRUE)
#'
#' # Information values of an item pool with multiple models
#' ip <- generate_ip(model = c("2PL", "3PL", "GPCM", "3PL", "GPCM"))
#' theta <- rnorm(sample(6:10, 1))
#' info(ip = ip, theta = theta[1])
#' info(ip = ip, theta = theta)
#' info(ip = ip, theta = theta, tif = TRUE)
#'
setMethod(f = "info", signature = c(ip = "Itempool"), .info_itempool)


###############################################################################@
############################# info (Testlet) ###################################
###############################################################################@
#' @export
#' @rdname info
#'
#' @examples
#' t1 <- generate_testlet(item_models = c("2PL", "3PL", "GRM", "3PL", "GRM"))
#' theta <- rnorm(sample(6:10, 1))
#' info(ip = t1, theta = theta[1])
#' info(ip = t1, theta = theta)
#' info(ip = t1, theta = theta, tif = TRUE)
setMethod(f = "info", signature = c(ip = "Testlet"), .info_itempool)


###############################################################################@
############################# info (REST) ######################################
###############################################################################@
#' @export
#' @rdname info
#'
setMethod(
  f = "info", signature = c(ip = "numMatDfListChar"),
  function(ip, theta, tif = FALSE, observed = FALSE, resp = NULL){
    if (inherits(ip, c("numeric", "data.frame", "matrix", "list"))) {
      return(info(ip = itempool(ip), theta = theta, tif = tif))
    } else {
      stop("Cannot convert object to an 'Item' or an 'Itempool' object. ",
           "Please provide a valid 'Item' or 'Itempool' object using either ",
           "'item()' or 'itempool()' function.")
    }
  }
)














#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% info_kl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


###############################################################################@
############################# info_kl (generic) ################################
###############################################################################@
#' Calculate Kullback-Leibler Information
#'
#' @description
#' \code{info_kl} returns the the value of Kullback-Leibler Information function.
#' The formulas are implemented from Chang and Ying (1996).
#'
#'
#' @param ip An \code{\link{Item-class}}, \code{\link{Itempool-class}} or a
#'   \code{\link{Testlet-class}} object.
#' @param trueTheta True theta parameter
#' @param theta Estimated theta parameter
#'
#' @return A vector (or matrix) consist of item or test Kullback-Leibler
#' information.
#'
#' @include Item-class.R
#' @include Itempool-class.R
#' @include Item-class-methods.R
#' @include Itempool-class-methods.R
#'
#' @author Emre Gonulates
#'
#' @keywords internal
#'
#' @references
#' Chang, H.-H., & Ying, Z. (1996). A Global Information Approach to
#' Computerized Adaptive Testing. Applied Psychological Measurement, 20(3),
#' 213-229. doi: 10.1177/014662169602000303
#'
#' #rdname info_kl
#' @noRd
#'
setGeneric("info_kl", function(ip, trueTheta, theta)
           {standardGeneric ("info_kl")})


###############################################################################@
############################# info_kl (Item) ###################################
###############################################################################@
#'
#' #rdname info_kl
#'
#' @noRd
#'
#' @keywords internal
#'
setMethod(
  f = "info_kl", signature = c(ip = "Item"),
  function(ip, trueTheta, theta){

    # Examples:
    # ip <- item(a = 3, b = 0, c = .1)
    # info_kl(ip = ip, trueTheta = 0, theta = 0)
    # info_kl(ip = ip, trueTheta = 1, theta = 1)
    # info_kl(ip = ip, trueTheta = c(-4, 4), theta = -4)
    #
    # # Reproduction of Figure 1 of Chang and Ying (1996).
    # persp(x = seq(-4,4,.1), y = seq(-4,4,.1),
    #       z = info_kl(ip = ip, trueTheta = seq(-4,4,.1), theta = seq(-4,4,.1)),
    #       xlab = "True Theta", ylab = "Estimated Theta",
    #       zlab = "Kullback-Leibler Information",
    #       phi = 45, theta = -45, ticktype = "detailed")
    #

    # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
    if (ip$model %in%
          names(Pmodels)[sapply(Pmodels, function(x) x$model_family == 'UIRT')])
    {
      if ((length(trueTheta) > 1) || (length(theta) > 1))
      {
        pTrue <- matrix(prob(ip = ip, theta = trueTheta),
                        nrow = length(trueTheta), ncol = length(theta))
        p <- matrix(prob(ip = ip, theta = theta),
                    nrow = length(trueTheta), ncol = length(theta), byrow = T)
      } else {
        pTrue <- prob(ip = ip, theta = trueTheta)
        p <- prob(ip = ip, theta = theta)
      }
      return(pTrue * log(pTrue / p) + (1-pTrue) * log((1-pTrue)/(1-p)))
    } else stop("This model is not implemented in this function yet.")
  }
)

