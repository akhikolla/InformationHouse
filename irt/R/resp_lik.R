
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% resp_lik %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############################################################################@
############################# resp_lik (generic) ###############################
###############################################################################@
#' Likelihood of a response string
#'
#' @description
#' \code{resp_lik} returns the likelihood of a response string
#' for given items and ability.
#'
#' @param resp A vector of item responses.
#' @param ip An \code{\link{Item-class}}, \code{\link{Itempool-class}} or a
#'   \code{\link{Testlet-class}} object.
#' @param theta An vector containing ability parameters.
#'
#' @return A matrix of likelihood(s)
#'
#' @include Item-class.R
#' @include Itempool-class.R
#' @include Item-class-methods.R
#' @include Itempool-class-methods.R
#'
#' @author Emre Gonulates
#'
#' @rdname resp_lik
#'
setGeneric("resp_lik", function(ip, resp, theta)
{standardGeneric ("resp_lik")})


###############################################################################@
############################# resp_lik (Item) ##################################
###############################################################################@
#' @export
#'
#' @rdname resp_lik
#'
#' @examples
#' item <- generate_item(model = "3PL")
#' theta <- rnorm(6)
#' resp <- sim_resp(ip = item, theta = theta, prop_missing = .1)
#' resp_lik(ip = item, resp = resp, theta = theta)
#'
#' item <- generate_item(model = "GRM")
#' resp <- sim_resp(ip = item, theta = theta, prop_missing = .1)
#' resp_lik(ip = item, resp = resp, theta = theta)
setMethod(
  f = "resp_lik", signature = c(ip = "Item"),
  function(ip, resp, theta){
    # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
    if (ip$model %in% c(names(Pmodels)[sapply(
      Pmodels, function(x) x$model_family == 'UIRT')], "GRM", "GPCM")) {
      if (inherits(resp, 'integer')) {
        if (all(is.na(resp))) return(resp)
        if (!all(is.na(resp) | resp %in% 0L:length(ip$b)))
          stop(paste0("\nInvalid response. Response should be ",
                      "either: ", paste0(0L:length(ip$b), collapse = ", "),
                      " or NA (i.e. missing)."))
        if (length(theta) == length(resp)) {
          return(resp_lik_item_cpp(resp = resp, theta = theta, item = ip))
        } else
          stop(paste0("Invalid arguments. Number of subjects (theta) ",
                      "should be equal to the number of responses."))
      } else if (inherits(resp, c("numeric", "logical"))) {
        return(resp_lik(ip = ip, resp = as.integer(resp), theta = theta))
      } else if (inherits(resp, c("matrix"))) {
        if ((ncol(resp) == 1) && (nrow(resp) == length(theta)))
        {
          return(resp_lik(ip = ip, resp = as.integer(resp), theta = theta))
        } else
          stop(paste0("\nInvalid arguments. Either resp has more than ",
                      "1 column, or number of rows in resp matrix is ",
                      "not equal to the length of theta vector."))
      } else if (inherits(resp, c("data.frame"))) { # Deal with tibbles as well
        return(resp_lik(ip = ip, resp = as.matrix(resp), theta = theta))
      } else
        stop(paste0("Invalid response. Response cannot be a ", class(resp)[1],
                    " object"))
    } else stop("\nThis model is not implemented in this function.")
  }
)


###############################################################################@
############################# resp_lik (Itempool) #############################
###############################################################################@
#' @export
#' @rdname resp_lik
#'
#' @examples
#' ip <- generate_ip(model = "3PL")
#' theta <- rnorm(6)
#' resp <- sim_resp(ip = ip, theta = theta, prop_missing = .1)
#' resp_lik(ip = ip, resp = resp, theta = theta)
#'
#' ip <- generate_ip(model = "GRM")
#' resp <- sim_resp(ip = ip, theta = theta, prop_missing = .1)
#' resp_lik(ip = ip, resp = resp, theta = theta)
#' 
setMethod(
  f = "resp_lik", signature = c(ip = "Itempool"),
  function(ip, resp, theta){
    # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM",
    #                                       "GRM", "GPCM"
    if (all(ip$model %in% c(names(Pmodels)[sapply(
      Pmodels, function(x) x$model_family == 'UIRT')], "GRM", "GPCM") |
      sapply(ip@item_list, is, "Testlet")))
    {
      noItem <- get_itempool_size(ip)["items"]
      noTheta <- length(theta)
      if (inherits(resp, 'integer')) {
        if (noItem == 1) {
          # This Part newly added FROM HERE
          if (noTheta == length(resp))
          {
            item <- ip@item_list[[1]]
            if (is(item, "Item")) {
              result <- resp_lik_item_cpp(resp = resp, theta = theta,
                                          item = item)
            } else if (is(item, "Testlet")) {
              result <- resp_lik_bare_testlet_cpp(resp = resp, theta = theta,
                                                  testlet = item)
            }
          } else
            stop(paste0("Invalid arguments. Number of subjects (theta) ",
                        "should be equal to the number of responses."))
        } else {
          if (all(is.na(resp))) return(NA)
          if ((noTheta == 1) && (noItem == length(resp))) {
            result <- resp_lik_itempool_cpp(resp = matrix(resp, nrow = 1),
                                             theta = theta, ip = ip)
            } else
              stop(paste0("Invalid arguments. Number of subjects (theta) ",
                          "should be equal to 1 and the number of ",
                          "responses should equal to the number of items."))
        }
      } else if (inherits(resp, c("numeric", "logical"))) {
        return(resp_lik(resp = as.integer(resp), ip = ip, theta = theta))
      } else if (inherits(resp, c("matrix"))) {
        if ((ncol(resp) == noItem) && (nrow(resp) == noTheta) &&
            is.numeric(resp)) {
          result <- resp_lik_itempool_cpp(resp = resp, theta = theta, ip = ip)
        } else
          stop(paste0("\nInvalid arguments. Number of columns of the resp ",
                      "matrix should be equal to the number of items.",
                      "Also number of rows of the resp matrix should be ",
                      "equal to the number of subjects (theta)."))
      } else if (inherits(resp, c("data.frame"))) {# Deal with tibbles as well
        result <- resp_lik(ip = ip, resp = sapply(resp, as.integer),
                           theta = theta)
      } else stop(paste0("Invalid response. Response cannot be a ",
                         class(resp)[1], " object"))
      if (!is.null(names(theta))) names(result) <- names(theta)
      return(result)
    } else stop(paste0("There is an item model in the item pool that is ",
                       "not implemented in this function yet."))
  }
)


###############################################################################@
############################# resp_lik (Testlet) ###############################
###############################################################################@
#' @export
#' @rdname resp_lik
#'
setMethod(
  f = "resp_lik", signature = c(ip = "Testlet"),
  function(ip, resp, theta){
    return(resp_lik(ip = itempool(ip), resp = resp, theta = theta))
  }
)

