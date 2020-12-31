setClassUnion(name = "listNULL", members = c("list","NULL"))
setClassUnion(name = "characterNULL", members = c("character","NULL"))
setClassUnion(name = "numericNULL", members = c("numeric","NULL"))
setClassUnion(name = "mNumeric", members = c("missing", "numeric"))
setClassUnion(name = "mCharacter", members = c("missing", "character"))
setClassUnion(name = "mLogical", members = c("missing", "logical"))
# Create a class union for: numeric, matrix, data.frame, list
setClassUnion(name = "numMatDfListChar",
              members = c("numeric", "character", "matrix", "data.frame",
                          "list"))

################################################################################
############################# Item class ######################################@
###########################################################################@####
#' An S4 class to represent an Item
#' @description
#' \code{Item} is a class to represent an item. An object in Item class should
#' have a \code{model} name and \code{parameters}.
#'
#' @slot id Item id. Default value is \code{NULL}.
#' @slot model The model that item \code{parameters} represents. Currently,
#'   following models are available:
#'     \describe{
#'       \item{\code{"Rasch"}}{
#'         Rasch Model.
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"b"}}{Item difficulty parameter.}
#'           }
#'
#'         Probability of correct response at ability estimate \eqn{\theta}:
#'
#'         \deqn{P(\theta) = \frac{e^{(\theta - b)}}{1+e^{(\theta - b)}}}
#'
#'         Model family: Unidimensional Item Response Theory (UIRT) Models
#'         }
#'       \item{\code{"1PL"}}{
#'         Unidimensional One-Parameter Logistic Model.
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"b"}}{Item difficulty parameter.}
#'           \item{\code{"D"}}{Scaling constant. Default value is \code{1}.}
#'           }
#'
#'         Probability of correct response at ability estimate \eqn{\theta}:
#'
#'         \deqn{P(\theta) = \frac{e^{D(\theta - b)}}{1+e^{D(\theta - b)}}}
#'
#'         Model family: Unidimensional Item Response Theory (UIRT) Models
#'         }
#'       \item{\code{"2PL"}}{
#'         Unidimensional Two-Parameter Logistic Model.
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"a"}}{Item discrimination parameter.}
#'           \item{\code{"b"}}{Item difficulty parameter.}
#'           \item{\code{"D"}}{Scaling constant. Default value is \code{1}.}
#'           }
#'
#'         Probability of correct response at ability estimate \eqn{\theta}:
#'
#'         \deqn{P(\theta) = \frac{e^{Da(\theta - b)}}{1+e^{Da(\theta - b)}}}
#'
#'         Model family: Unidimensional Item Response Theory (UIRT) Models
#'         }
#'       \item{\code{"3PL"}}{
#'         Unidimensional Three-Parameter Logistic Model.
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"a"}}{Item discrimination parameter.}
#'           \item{\code{"b"}}{Item difficulty parameter.}
#'           \item{\code{"c"}}{Pseudo-guessing parameter (lower asymptote).}
#'           \item{\code{"D"}}{Scaling constant. Default value is \code{1}.}
#'           }
#'
#'         Probability of correct response at ability estimate \eqn{\theta}:
#'
#'         \deqn{P(\theta) = c + (1-c) \frac{e^{Da(\theta - b)}}{1+e^{Da(\theta - b)}}}
#'
#'         Model family: Unidimensional Item Response Theory (UIRT) Models
#'         }
#'       \item{\code{"4PL"}}{
#'         Unidimensional Four-Parameter Logistic Model.
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"a"}}{Item discrimination parameter.}
#'           \item{\code{"b"}}{Item difficulty parameter.}
#'           \item{\code{"c"}}{Pseudo-guessing parameter (lower asymptote).}
#'           \item{\code{"d"}}{Upper asymptote parameter.}
#'           \item{\code{"D"}}{Scaling constant. Default value is \code{1}.}
#'           }
#'
#'         Probability of correct response at ability estimate \eqn{\theta}:
#'
#'         \deqn{P(\theta) = c + (d-c) \frac{e^{Da(\theta - b)}}{1+e^{Da(\theta - b)}}}
#'
#'         Model family: Unidimensional Item Response Theory (UIRT) Models
#'         }
#'       \item{\code{"GRM"}}{
#'         Graded Response Model
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"a"}}{Item discrimination parameter.}
#'           \item{\code{"b"}}{Item threshold parameters (a vector of values).
#'             Each value refers to the ability level for which the probability
#'             of responding at or above that category is equal to 0.5. }
#'           \item{\code{"D"}}{Scaling constant. Default value is \code{1}.}
#'           }
#'
#'         Probability of scoring at or above the category \eqn{k}:
#'
#'         \deqn{P^*_k(\theta) = \frac{e^{Da(\theta - b_k)}}{1+e^{Da(\theta - b_k)}}}
#'
#'         Probability of responding at category \eqn{k} where the possible
#'         scores are \eqn{0, \ldots, m}:
#'
#'         \deqn{P_0(\theta) = 1 - P^*_1(\theta)}
#'         \deqn{P_1(\theta) = P^*_1(\theta) - P^*_2(\theta)}
#'         \deqn{\cdots}
#'         \deqn{P_k(\theta) = P^*_{k}(\theta) - P^*_{k+1}(\theta)}
#'         \deqn{\cdots}
#'         \deqn{P_m(\theta) = P^*_{m}(\theta)}
#'
#'         Model family: Polytomous Item Response Theory (PIRT) Models
#'         }
#'       \item{\code{"GPCM"}}{
#'         Generalized Partial Credit Model
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"a"}}{Item discrimination parameter.}
#'           \item{\code{"b"}}{Item step difficulty parameters (a vector of
#'             values).}
#'           \item{\code{"D"}}{Scaling constant. Default value is \code{1}.}
#'           }
#'
#'         Probability of scoring at category \eqn{k}:
#'
#'         \deqn{P_k(\theta) = \frac{exp[\sum_{v = 0}^{k} Da(\theta - b_v)]}
#'           {\sum_{c = 0}^{m-1}exp[\sum_{v = 0}^{c}Da(\theta - b_v)]}}
#'
#'         Model family: Polytomous Item Response Theory (PIRT) Models
#'         }
#'       \item{\code{"PCM"}}{
#'         Partial Credit Model (Masters, 1982)
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"b"}}{Item step difficulty parameters (a vector of
#'             values).}
#'           }
#'
#'         Probability of scoring at category \eqn{k}:
#'
#'         \deqn{P_k(\theta) = \frac{exp[\sum_{v = 0}^{k} (\theta - b_v)]}{\sum_{c = 0}^{m-1}exp[\sum_{v = 0}^{c}(\theta - b_v)]}}
#'
#'         Model family: Polytomous Item Response Theory (PIRT) Models
#'         }
#'       \item{\code{"GPCM2"}}{
#'         An alternative parametrization of Generalized Partial Credit Model
#'         \code{"GPCM"} where \eqn{b_k = b - d_k}. See Muraki (1997),
#'         Equation 15 on page 164.
#'
#'         Required parameters:
#'
#'         \describe{
#'           \item{\code{"a"}}{Item discrimination parameter.}
#'           \item{\code{"b"}}{Location parameter.}
#'           \item{\code{"d"}}{A vector of threshold parameters.}
#'           \item{\code{"D"}}{Scaling constant. Default value is \code{1}.}
#'           }
#'
#'         Probability of scoring at category \eqn{k}:
#'
#'         \deqn{P_k(\theta) = \frac{exp[\sum_{v = 0}^{k} Da(\theta - b + d_v)]}{\sum_{c = 0}^{m-1}exp[\sum_{v = 0}^{c}Da(\theta - b + d_v)]}}
#'
#'         Model family: Polytomous Item Response Theory (PIRT) Models
#'         }
#'     }
#'
#'   A model must be specified for the construction of an \code{Item} object.
#'
#' @slot parameters A list containing numeric vectors that represent item
#'   parameters. Depending on the model these can change.
#' @slot se_parameters Standard error of the item parameters. This should be
#'   a list of standard error values. For example, for "2PL", if the parameters
#'   are \code{list(a = 1.2, b = -0.22)}, the standard error values of
#'   parameters can be either \code{NULL} (which is the default value) or
#'   \code{list(a = 0.24, b = 0.42)}. None of the standard error values can
#'   be smaller than 0. Individual SE values can be \code{NA}. For example,
#'   \code{list(a = 0.24, b = NA)} is acceptable, whereas
#'   \code{list(a = 0.24, b = NULL)} is not acceptable.
#'
#'   For models like polytomous items, the SE values should match the parameter
#'   values in length. For example, if the parameter values of a \code{"GPCM"}
#'   is \code{parameters = list(a = 1.4, b = c(-1, 0.42, 2.1), D = 1.7)}, then
#'   the SE values should be like
#'   \code{se_parameters = list(a = .2, b = c(.32, 0.34, .3))}. Since the
#'   scaling parameter \code{D} is constant, it does not have a standard error.
#'
#' @slot content Content information for the Item object.
#' @slot misc This slot is a list where one can put any information about
#'  the Item object. For example, one can enter the id's of the enemies of the current
#'  Item as \code{misc = list(enemies = c("i1", i2))}. Or, one can enter
#'  Sympson-Hetter exposure control parameter K:
#'  \code{misc = list(sympson_hetter_k = .75)}.
#'
#' @export
#'
#' @references
#'
#'   Masters, G. N. (1982). A Rasch model for partial credit scoring.
#'   \emph{Psychometrika}, 47, 149–174.
#'
#'   Muraki, E. (1992). A generalized partial credit model:
#'   Application of an EM algorithm. \emph{Applied Psychological Measurement},
#'   16, 159–176.
#'
#'
setClass(Class = "Item",
         slots = c(id = "characterNULL",
                   model = "characterNULL",
                   parameters = "listNULL",
                   se_parameters = "listNULL",
                   content = "characterNULL",
                   misc = "listNULL")
         )

# The default scaling parameter value
default_D_scaling <-  1 # 1.702
default_D_scaling_mirt <- 1 # 1.702

# List of currently implemented models for Item class
# Pmodels : Psychometric Models
# Get all Unidimensioal IRT (UIRT) models
# names(Pmodels)[sapply(Pmodels, function(x) x$model_family == 'UIRT')]
Pmodels <- list(
  'Rasch' = list(parameters = list(
    b = list(name = 'b', size = 1, se = TRUE,
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0)),
    model_family = "UIRT",
    verbose_name = "Rasch Model"),
  '1PL' = list(parameters = list(
    b = list(name = 'b', size = 1, se = TRUE,
             description = "Difficulty Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling)),
    model_family = "UIRT",
    verbose_name = "Unidimensional One-Parameter Logistic Model"),
  '2PL' = list(parameters = list(
    a = list(name = 'a', size = 1, se = TRUE,
             description = "Discrimination Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 1),
    b = list(name = 'b', size = 1, se = TRUE,
             description = "Difficulty Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling)),
    model_family = "UIRT",
    verbose_name = "Unidimensional Two-Parameter Logistic Model"),
  '3PL' = list(parameters = list(
    a = list(name = 'a', size = 1, se = TRUE,
             description = "Discrimination Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 1),
    b = list(name = 'b', size = 1, se = TRUE,
             description = "Difficulty Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0),
    c = list(name = 'c', size = 1, se = TRUE, min = 0, max = 1,
             description = "Pseudo-Guessing Parameter",
             default_value = 0),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling)),
    model_family = "UIRT",
    verbose_name = "Unidimensional Three-Parameter Logistic Model"),
  '4PL' = list(parameters = list(
    a = list(name = 'a', size = 1, se = TRUE,
             description = "Discrimination Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 1),
    b = list(name = 'b', size = 1, se = TRUE,
             description = "Difficulty Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0),
    c = list(name = 'c', size = 1, se = TRUE, min = 0, max = 1,
             description = "Pseudo-Guessing Parameter",
             default_value = 0),
    d = list(name = 'd', size = 1, se = TRUE, min = 0, max = 1,
             description = "Upper Asymptote Parameter",
             default_value = 1),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling)),
    model_family = "UIRT",
    verbose_name = "Unidimensional Four-Parameter Logistic Model"),
  'M1PL' = list(parameters = list(
    d = list(name = 'd', size = 1, se = TRUE,
             description = "Intercept Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling_mirt)),
    model_family = "MIRT",
                verbose_name = "Multidimensional One-Parameter Logistic Model"),
  'M2PL' = list(parameters = list(
    a = list(name = 'a', size = 100, se = TRUE,
             description = "Slope Parameters",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 1),
    d = list(name = 'd', size = 1, se = TRUE,
             description = "Intercept Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling_mirt)),
    model_family = "MIRT",
    verbose_name = "Multidimensional Two-Parameter Logistic Model"),
  'M3PL' = list(parameters = list(
    a = list(name = 'a', size = 100, se = TRUE,
             description = "Slope Parameters",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 1),
    c = list(name = 'c', size = 1, se = TRUE, min = 0, max = 1,
             description = "Pseudo-Guessing Parameter",
             default_value = 0),
    d = list(name = 'd', size = 1, se = TRUE,
             description = "Intercept Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling_mirt)),
    model_family = "MIRT",
    verbose_name = "Multidimensional Three-Parameter Logistic Model"
                ),
  'GRM' = list(parameters = list(
    a = list(name = 'a', size = 1, se = TRUE,
             description = "Discrimination Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 1),
    b = list(name = 'b', size = 100, se = TRUE,
             description = "Threshold Parameters",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = c(-1, 0, 1)),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling)),
    model_family = "PIRT",
    verbose_name = "Graded Response Model"),
  'PCM' = list(parameters = list(
    b = list(name = 'b', size = 100, se = TRUE,
             description = "Step Difficulty Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = c(-1, 0, 1))),
    model_family = "PIRT",
    verbose_name = "Partial Credit Model"),
  'GPCM' = list(parameters = list(
    a = list(name = 'a', size = 1, se = TRUE,
             description = "Discrimination Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 1),
    b = list(name = 'b', size = 100, se = TRUE,
             description = "Step Difficulty Parameters",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = c(-1, 0, 1)),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling)),
    model_family = "PIRT",
    verbose_name = "Generalized Partial Credit Model"),
  'GPCM2' = list(parameters = list(
    a = list(name = 'a', size = 1, se = TRUE,
             description = "Discrimination Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 1),
    b = list(name = 'b', size = 1, se = TRUE,
             description = "Overall Location Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = 0),
    d = list(name = 'd', size = 100, se = TRUE,
             description = "Threshold Parameters",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = c(-1, 0, 1)),
    D = list(name = 'D', size = 1, se = FALSE,
             description = "Scaling Parameter",
             min = -.Machine$double.xmax, max = .Machine$double.xmax,
             default_value = default_D_scaling)),
    model_family = "PIRT",
    verbose_name = "Reparametrized Generalized Partial Credit Model")
  )



################################################################################
############################# initialize (Item) ###############################@
###########################################################################@####
#' @noRd
#' @title This function initializes the \code{Item} object.
#'
#' @importFrom methods callNextMethod
#'
setMethod("initialize", "Item",
          function(.Object, id = NULL, model = "2PL",
                   parameters = list(a = 1, b = 0, D = default_D_scaling),
                   se_parameters = NULL, content = NULL, misc = NULL, ...) {
  .Object <- callNextMethod(.Object, ...)
  .Object@id <- id
  if (is.character(model) && model %in% names(Pmodels))
    .Object@parameters <- parameters else
      stop(paste0("Invalid model value. Model name should be one of the ",
                  'following:\n"', paste0(names(Pmodels), collapse = '", "'),
                  '"'))
  .Object@model <- model
  if (is.list(parameters)) .Object@parameters <- parameters else
    stop("Invalid parameters. Parameters should be a 'list'.")
  .Object@se_parameters <- se_parameters
  .Object@content <- content
  .Object@misc <- misc
  # Check validity of the object
  validObject(.Object)
  # The parameters is a list and the user can give the parameter names in any
  # order. Here, this function fixes the order to the order seen in
  # Pmodels[['Rasch']]$parameters.
  .Object@parameters <- parameters[names(Pmodels[[model]]$parameters)]
  .Object
})



################################################################################
############################# setValidity (Item) ##############################@
###########################################################################@####
# Here are some rules
# Each Item object should:
# * Have model and parameters slot
# * Have specified number of parameters. No more, no less.
# * Have proper names for parameters: all names should be listed in
#   names(Pmodels[[model_name]]$parameters) and should be unique.
# * Length of each sub-parameter should be 1 for IRTDichotomous models.
# * Each sub-parameter should be "numeric"
# * not be NULL or NA
# * not include NA in parameters slot
#
# TO-DO
# * se_parameters should also be named accordingly. And relevant checks should
#   be done to conform with parameter values.
#' @noRd
#' @title This function sets some validity rules for \code{Item} object.
#'
#' @name Item-class::setValidity()
#'
#' @param Class The class of the object.
#'
#'
setValidity(
  Class = "Item",
  function(object)
  {
    ############# LOCAL FUNCTIONS #############################################@
    # This function finds the parName (parameter names) within the
    # pars (parameters)
    # @description
    #
    # @param pars Parameter vector. Elements can have names or not.
    # @param parName A chacacter vector with length 1. Parameter name to
    # be searched within the names of \code{pars}.
    # @param irtModel whether the pars are the IRT model or not.
    # Specifically 1-4 parameter IRT models. For this function to work when
    # irtModel is TRUE, parName should be either "a", "b", "c" or "d".
    #
    # @return the value of the parameter searched.
    #
    # @examples
    # findParameter(c(a=1, b=2,c=.1, d=.9), "d")
    # findParameter(c(1, 2, .1, .9), "c")
    # findParameter(list(a = 2, b = .22), "a")
    # findParameter(list(a = 2, b = .22, c = .2), c("a", "c"))
    # findParameter(list(a = c(1.2, 2), b = .22, c = .2), "a")
    findParameter <- function(pars, parName, irtModel = TRUE)
    {
      if (any(names(pars) %in% parName) &&
          (sum(names(pars) %in% parName) == 1))
      {
        result <- pars[names(pars) %in% parName]
      } else if (irtModel)
        result <- pars[which(letters %in% parName)]
      return(unlist(result))
    }

    checkCParameter <- function(pars)
    {
      # Find c parameter from the list of pars
      # Examples:
      # checkCParameter(object@parameters)
      cPar <- findParameter(pars, "c", irtModel = TRUE)
      if ((cPar < 0) || (cPar > 1))
        stop(paste0("Invalid 'c' parameter. 'c' parameter in IRT ",
                    "parameterization cannot be smaller than 0 or larger ",
                    "than 1. Please check your 'c' parameter"), call. = FALSE)
    }
    ############## End of LOCAL FUNCTIONS #####################################@

    # ----------------------- Check model ------------------------------- #
    # Check the model, currently only irt1PM, irt2PM,
    # irt3PM, irt4PM, mirt1PM, mirt2PM and mirt3PM is used
    if (is.null(object@model) || (length(object@model) != 1) ||
          !(object@model %in% names(Pmodels)))
      stop(paste0("Invalid model. Item model should be specified ",
                  "correctly. It can be either: ",
                  paste0(names(Pmodels), collapse = ", ")), call. = FALSE)
    # This will signify the model name in Pmodels
    model_name <- object@model

    # ----------------------- Check parameters ------------------------------- #
    # Object parameters cannot be NULL or NA
    if (is.null(object@parameters) || any(is.na(object@parameters)))
      stop("Invalid parameter. Item parameters cannot be NULL or NA.",
           call. = FALSE)

    # All parameters should be numeric
    if (!all(sapply(object@parameters, "class") %in% c("integer", "numeric")))
      stop(paste0("Invalid parameters. All parameters should be numeric."),
           call. = FALSE)

    # Check for proper naming of parameters. Parameter names should be unique
    # and all should correspond one-to-one with Pmodels[[model_name]]$parameters
    if (is.null(parNames <- names(object@parameters))) {
      stop(paste0("Invalid parameter names. Parameter names of Item class ",
                  "cannot be NULL. Please give relevant names."), call. = FALSE)
    } else if ((
      length(parNames) !=  length(Pmodels[[model_name]]$parameters)) ||
      length(unique(parNames)) != length(Pmodels[[model_name]]$parameters)
      ) {
      stop(paste0("Invalid parameter names. Parameter names of Item class ",
                  "should be unique and complete. Please give relevant names."),
           call. = FALSE)
    } else if (!all(parNames %in% names(Pmodels[[model_name]]$parameters)))
      stop(paste0("Invalid parameter names. Parameter names for ",
                  model_name," model should be ",
                  paste0(names(Pmodels[[model_name]]$parameters),
                         collapse = ", "),
                  ". Please give relevant names."), call. = FALSE)

    # Number of parameters should be as specified in Pmodels parameters
    if (length(object@parameters) != length(Pmodels[[model_name]]$parameters))
      stop(paste0("Invalid parameters. Number of parameters for ",
                  model_name," model should be ",
                  length(Pmodels[[model_name]]$parameters), "."), call. = FALSE)

    # Size of the model parameters should be matching to the sizes designated
    # in Pmodels.
    # Also the magnitudes of parameter values should be within min-max values
    # designated in Pmodels
    for (p in Pmodels[[model_name]]$parameters) {
      temp_par <- object@parameters[[p$name]]
      if (length(temp_par) < 1 || length(temp_par) > p$size)
        stop(paste0("Invalid parameters. For \"", model_name,"\" model the ",
                    "size of the \"", p$name, "\" should be ",
                    ifelse(p$size == 1, 1, paste0(p$size, " or less")), "."),
             call. = FALSE)
      if (any(temp_par < p$min) || any(temp_par > p$max))
        stop(paste0("Invalid parameters. The values of \"", p$name, "\" ",
                    "parameters should be between ", prettyNum(p$min), " and ",
                    prettyNum(p$max), "."), call. = FALSE)
    }
    if (Pmodels[[model_name]]$model_family == "UIRT") {
      # Length of all sub-parameters should be 1 for
      if (any(sapply(object@parameters, FUN = "length") != 1))
        stop(paste0("Invalid parameters. For Unidimensional IRT models, the ",
                    "length of all parameters should be 1."), call. = FALSE)
      switch(object@model,
             "3PL" = {
               checkCParameter(object@parameters)
             },
             "4PL" = {
               checkCParameter(object@parameters)
               checkDParameter <- function(pars)
               {
                 cPar <- findParameter(pars, "c", irtModel = TRUE)
                 dPar <- findParameter(pars, "d", irtModel = TRUE)
                 if ((dPar < 0) || (dPar > 1) || (dPar < cPar))
                   stop(paste0("Invalid 'd' parameter. 'd' parameter in IRT ",
                               "parameterization cannot be smaller than 0, ",
                               "larger than 1 and smaller than 'c' parameter. ",
                               "Please check your 'd' parameter."),
                        call. = FALSE)
               }
               checkDParameter(object@parameters)
             }
      )
    } else if (Pmodels[[model_name]]$model_family == "MIRT") {
      # d parameter should have length 1.
      if (length(object@parameters$d) != 1)
        stop(paste0("Invalid parameters. The length of 'd' parameter should ",
                    "be 1."), call. = FALSE)
      # Check c parameter
      if (object@model == "M3PL")
      {
        # c parameter should have length 1.
        if (length(object@parameters$c) != 1)
          stop(paste0("Invalid parameters. The length of 'c' parameter should ",
                      "be 1."), call. = FALSE)
        # c parameter should be between 0 and 1.
        checkCParameter(object@parameters)
      }
    } else if (Pmodels[[model_name]]$model_family == "PIRT") {
      # item location parameter "b" cannot be a vector of length more than 1
      if (object@model == "GPCM2") {
        if (length(object@parameters$b) != 1)
          stop("In 'GPCM2' model, the item location parameter 'b' should be ",
               "a single value (i.e. length(b) should be 1).", call. = FALSE)
      }

    }

    # ----------------------- Check se_parameters ---------------------------- #
    if (!is.null(object@se_parameters)) {
      se_par_names <- unlist(sapply(Pmodels[[model_name]]$parameters,
                                    function(x) if (x$se) x$name))

      if (is.null(names(object@se_parameters)) ||
          !all(names(object@se_parameters) %in% se_par_names)||
          !all(se_par_names %in% names(object@se_parameters))
          )
        stop(paste0("Invalid 'se_parameters' values.\n'se_parameters' should ",
                    "be a list with elements named: ",
                    paste0("'", se_par_names, "'", collapse = ", "), "."),
             call. = FALSE)
      if (length(object@se_parameters) != length(se_par_names))
        stop(paste0("Invalid 'se_parameters' values. 'se_parameters' length ",
                    "for ", model_name," model should be ",
                    length(se_par_names), "."), call. = FALSE)

      # Individual se_parameters cannot be NULL
      if (all(any(sapply(object@se_parameters, is.null))))
        stop(paste0("Invalid 'se_parameters' values. Individual elements of ",
                    "'se_parameters' cannot be NULL."),
             call. = FALSE)

      # The length of each se_parameter should be the same as parameters, i.e.
      # for example for "GPCM", if there are 3 threshold parameters, there
      # should be three standard errors for each threshold parameter.
      if (!all(sapply(object@parameters[se_par_names], length)[se_par_names] ==
               sapply(object@se_parameters, length)[se_par_names]))
        stop(paste0("Invalid 'se_parameters' values. All of the elements of ",
                    "'se_parameters' should have the same length as the ",
                    "corresponding element in 'parameters'."),
             call. = FALSE)


      # All SE parameters should be numeric if it is not NULL or NA
      if (!all(sapply(object@se_parameters, FUN = function(x)
        all(is.na(x)) | (class(x) == "numeric"))))
        stop(paste0("Invalid 'se_parameters' values. All standard error values",
                    " of item parameters ('se_parameters') should be numeric."),
             call. = FALSE)
      # All se_parameters should be larger than 0 or they are NA
      if (any(!is.na(unlist(object@se_parameters)) &
              unlist(object@se_parameters) < 0))
        stop(paste0("Invalid 'se_parameters' values. Standard error values of ",
                    "item parameters ('se_parameters') cannot be smaller ",
                    "than 0."), call. = FALSE)
    }

    # ----------------------- Check id --------------------------------------- #
    # Check the length of id vector
    if (!is.null(object@id) && (length(object@id) != 1))
      stop(paste0("Invalid item id. Item id should have length 1, ",
                  "or be NULL."), call. = FALSE)
  })




# Per advice from:
# http://r-pkgs.had.co.nz/src.html#c-best-practices
.onUnload <- function (libpath) {
  library.dynam.unload("irt", libpath)
}
