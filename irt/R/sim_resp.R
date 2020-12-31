
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% sim_resp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


###############################################################################@
############################# sim_resp (generic) ###############################
###############################################################################@
#' Generate responses for a given model
#'
#' @description
#' \code{sim_resp} Generate dichotomous (0 or 1) or polytomous responses for
#' given ability and item parameter.
#'
#' @param ip An \code{\link{Item-class}}, \code{\link{Itempool-class}},
#'   \code{\link{Testlet-class}} object containing the item parameters.
#' @param theta An object containing the subject ability parameters.
#' @param prop_missing Proportion of responses that should be missing. Default
#'   value is \code{0}. This argument is valid for only
#'   \code{\link{Itempool-class}} and \code{\link{Testlet-class}} objects.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#' @include Itempool-class-methods.R
#'
#' @return A vector of responses.
#'
#' @author Emre Gonulates
#'
setGeneric("sim_resp", function(ip, theta, prop_missing = 0)
{standardGeneric ("sim_resp")})


###############################################################################@
############################# sim_resp (Item) ##################################
###############################################################################@
#' @export
#'
#' @rdname sim_resp
#'
#' @examples
#' ## Simulate Responses for an Item object ##
#' item <- generate_item(model = "3PL")
#' sim_resp(ip = item, theta = rnorm(1))
#'
#' item <- generate_item(model = "GPCM")
#' sim_resp(ip = item, theta = rnorm(1))
#'
#'
#' item <- generate_item(model = "GRM")
#' sim_resp(ip = item, theta = rnorm(1))
#'
setMethod(
  f = "sim_resp", signature = c(ip = "Item"),
  function(ip, theta, prop_missing = 0){
    # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
    if (ip$model %in%
          names(Pmodels)[sapply(Pmodels, function(x) x$model_family == 'UIRT')])
    {
      u <- runif(length(theta))
      P <- prob(ip = ip, theta = theta)
      return(stats::setNames(as.integer(u<P), names(theta)))
    } else if (ip$model %in% c("PCM", "GRM", "GPCM", "GPCM2")) {
      # This method is based on De Ayala (1994) The Influence of
      # Multidimensionality on the, p. 158-159: "The generation of an examinee’s
      # polytomous response string was accomplished by calculating the
      # probability of responding to each item alternative according to the
      # MGRM; the scaling factor D was set to 1.0. Based on the probability for
      # each alternative, cumulative probabilities were obtained for each
      # alternative. A random error component was incorporated into each
      # response by selecting a random number from a uniform distribution [0, 1]
      # and comparing it to the cumulative probabilities. The ordinal position
      # of the first cumulative probability that was greater than the random
      # number was taken as the examinee’s response to the item."
      u <- runif(length(theta))
      P <- prob(ip = ip, theta = theta)
      cP <- t(apply(P, MARGIN = 1, cumsum))
      # Subtract 1 so that the first category is 0.
      return(stats::setNames(apply(cP > u, 1, function(x) which(x)[1]) - 1L,
                             names(theta)))
    } else stop("This model is not implemented in this 'sim_resp' function.")
  }
)


###############################################################################@
############################# sim_resp (Testlet) ###############################
###############################################################################@
#' @export
#'
#' @rdname sim_resp
#'
#' @examples
#' ## Simulate Responses for a Testlet object ##
#' # Create a testlet
#' testlet <- testlet(c(item(b = 1), item(a = .8, b = 3.1),
#'                    item(b = -1:1, model = "PCM")))
#' sim_resp(ip = testlet, theta = rnorm(1))
setMethod(
  f = "sim_resp", signature = c(ip = "Testlet"),
  function(ip, theta, prop_missing = 0){
    ip <- ip@item_list
    return(sim_resp(ip, theta, prop_missing = prop_missing))
  }
)


###############################################################################@
############################# sim_resp (Itempool) #############################
###############################################################################@
#' @export
#'
#' @rdname sim_resp
#'
#' @examples
#' ## Simulate Responses for an Itempool object ##
#' # Create 3PL IRT item parameters
#' ip <- itempool(a = rlnorm(10, 0, 0.3), b = rnorm(10), c = runif(10, 0, .3))
#' # Simulate responses for one theta:
#' sim_resp(ip = ip, theta = rnorm(1))
#' # Simulate responses for eight thetas:
#' sim_resp(ip = ip, theta = rnorm(8))
#'
#' # Create Graded Response Model Parameters
#' ip <- generate_ip(n = 5, model = "GRM", n_categories = c(3, 4, 8, 5, 4))
#' # Simulate responses for one theta:
#' sim_resp(ip = ip, theta = rnorm(1))
#' # Simulate responses for 5 thetas:
#' sim_resp(ip = ip, theta = rnorm(5))
#' # Set 10% of the item responses as missing
#' sim_resp(ip = ip, theta = rnorm(5), prop_missing = .1)
setMethod(
  f = "sim_resp", signature = c(ip = "Itempool"),
  function(ip, theta, prop_missing = 0){
    # If the Model is one of the following: "irt1PM" "irt2PM" "irt3PM" "irt4PM"
    if (all(ip$model %in% c(names(Pmodels)[sapply(Pmodels, function(x)
      x$model_family == 'UIRT')], "GRM", "GPCM", "GPCM2", "PCM", "BTM")))
    {
      result <- sapply(ip@item_list, FUN =
                         function(x) sim_resp(ip = x, theta = theta))
      col_names <- ip$id
      if (!is.matrix(result)) {
        if (is.list(result)) { # If ip has testlets,then result should be a list
          col_names <- ip$resp_id
          result <- unlist(result)
          # return(matrix(unlist(result), nrow = length(theta),
          #               dimnames = list(names(theta), ip$resp_id)))
        }
        result <- matrix(result, nrow = length(theta))
        # result <- matrix(result, nrow = length(theta),
        #                  ncol = length(ip@item_list))
      }

      result[sample(1:length(result), round(length(result)*prop_missing))] <- NA
      # Set column and row names or resp matrix.
      if (is.null(names(theta))) {
        rownames(result) <- paste0("S", 1:length(theta))
        } else rownames(result) <- names(theta)
      colnames(result) <- col_names
      return(result)
    } else stop("This model is not implemented in 'sim_resp' function.")
  }
)


###############################################################################@
############################# sim_resp (REST) ##################################
###############################################################################@
#' @export
#'
#' @rdname sim_resp
#'
setMethod(
  f = "sim_resp", signature = c(ip = "numMatDfListChar"),
  function(ip, theta){
    if (inherits(ip, "numeric")) {
      tryCatch({
        return(sim_resp(ip = itempool(ip), theta = theta))
      }, error = function(e) {
        stop("Cannot convert object to an 'Item' object. Please provide a ",
             "valid object using 'item()' function. \nThe reason for ",
             "conversion failure: ", e)
      })
    } else if (inherits(ip, c("matrix", "data.frame", "list"))) {
      tryCatch({
        return(sim_resp(ip = itempool(ip), theta = theta))
      }, error = function(e) {
        stop("Cannot convert object to an 'Itempool' object. Please ",
             "provide a valid object using 'itempool()' function. \nThe ",
             "reason for conversion failure: \n", e)
      })
    } else
      stop("Cannot convert object to an 'Item' or an 'Itempool' object. ",
           "Please provide a valid 'Item' or 'Itempool' object using either ",
           "'item()' or 'itempool()' function.")
  }
)

