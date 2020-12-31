setClassUnion(name = "listNULL", members = c("list","NULL"))
setClassUnion(name = "characterNULL", members = c("character","NULL"))

################################################################################
############################# Itempool class #################################@
###########################################################################@####
#' An S4 class to represent an Itempool
#' @description
#' \code{\link{Itempool-class}} is a class to represent an item pool. This class is composed
#' of the collection of 'Item' class objects.
#'
#' @slot item_list The list of items that are 'Item' class
#' @slot misc A list of additional parameters for the item pool. For example,
#'   one can put the calibration date of the item pool as
#'   \code{misc = list(calibration_date = as.Date("2020-01-17"))}.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#'
#' @export
#'
#' @author Emre Gonulates
#'
setClass(
  Class="Itempool",
  slots = c(item_list="list", misc = "listNULL")
)


################################################################################
############################# initialize (Itempool) ##########################@
###########################################################################@####
#' @noRd
#' @title This function initializes an \code{Itempool} object.
#'
#' @importFrom methods callNextMethod new
#'
#' @author Emre Gonulates
#'
setMethod("initialize", "Itempool",
          function(.Object, item_list = list(new(Class = "Item")),
                   misc = NULL, ...) {
  .Object <- callNextMethod(.Object, ...)
  .Object@item_list <- item_list
  .Object@misc <- misc
  # Check validity of the object
  validObject(.Object)
  .Object
})



################################################################################
############################# setValidity(Itempool) ##########################@
###########################################################################@####
setValidity(
  Class = "Itempool",
  function(object)
  {
    # Itempool cannot be empty.
    if (length(object@item_list) == 0)
      stop("Invalid elements. Item pool cannot be empty. ")

    # All of the elements of the list should be "Item" or "Testlet" class.
    if (!all(sapply(object@item_list, FUN = function(x) is(x, "Item") | is(x, "Testlet")) ))
      stop(paste0("Invalid 'Itempool' elements. All of the elements of ",
                  "'Itempool' class should be either 'Item' or 'Testlet' ",
                  "class."))

    # Get first level Item id's, standalone items and Testlet items.
    item_ids <- unlist(lapply(object@item_list, FUN = function(x) x@id))
    # Make sure that the item_list has the names equal to the item_ids
    if (is.null(names(object@item_list)) ||
        !all(names(object@item_list) == item_ids))
      stop(paste0("The names of the 'item_list' elements should be the same ",
                  "as the id's of the 'Item' or 'Testlet' objects."))
    # Add the item id's of the testlet
    for (testlet in object@item_list[sapply(object@item_list, is, "Testlet")])
      item_ids <- c(item_ids, testlet@item_list$id)
    # id's of all elements should be unique and they cannot be NULL.
    if (is.null(item_ids) || any(duplicated(item_ids)))
      stop(paste0("Invalid id's. Each Item object in the item pool should have a ",
                  "unique id. Items in the testlets should not have
                  duplicate id's with standalone items."))
  })




################################################################################
############################# Testlet class ###################################@
###########################################################################@####
#' An S4 class to represent a Testlet
#' @description
#' \code{Testlet} is a class to represent an a collection of items. Items that
#' are connected by a common stimulus (for example a reading passage, a graph,
#' etc.) can form a testlet. An object in \code{Testlet} class should
#' have a \code{model} name and \code{item_list} which is an \code{Itempool}.
#' object. In fact, a \code{Testlet} object is very similar to an
#' \code{\link{Itempool-class}} object, except, it has a designated model and optional
#' parameters
#'
#' @slot id Testlet id. Default value is \code{NULL}.
#' @slot item_list A list of \code{Item} objects.
#' @slot model The model that testlet \code{parameters} represents. Currently
#' model can be:
#' BTM (Basic Testlet Model, this is default testlet model where no
#'      parameters necessary and testlet simply connects items),
#' RTM (Rasch Testlet Model),
#' BF (Bifactor Model) (Not implemented yet),
#' 2PTM (Two-parameter testlet model),
#' 3PTM (three-parameter testlet model).
#' A model must be specified for the construction of an \code{tetlet} object.
#' @slot parameters A list containing numeric vectors that represent testlet
#' parameters. Depending on the model these parameters can change.
#' @slot se_parameters Standard error of testlet parameters.
#' @slot content Content information for testlet.
#' @slot misc A list of additional parameters for the testlet.
#'
#' @export
#'
#' @author Emre Gonulates
#'
setClass(Class = "Testlet",
         slots = c(id = "characterNULL",
                   item_list="Itempool",
                   model = "character",
                   parameters = "listNULL",
                   se_parameters = "listNULL",
                   content = "characterNULL",
                   misc = "listNULL"))


# List of currently implemented models for Testlet class
# tmodels : Testlet Models
tmodels <- list(
  'BTM' = list(parameters= NULL,
               se_parameters = NULL,
               verbose_name = "Basic Testlet Model"),
  'RTM' = list(parameters= c('mean', 'sd'),
               se_parameters = c('mean', 'sd'),
               verbose_name = "Rasch Testlet Model")
  )


################################################################################
############################# initialize (Testlet) ############################@
###########################################################################@####
#' @noRd
#' @title This function initializes a \code{Testlet} object.
#'
#' @importFrom methods callNextMethod
#'
#' @author Emre Gonulates
#'
setMethod("initialize", "Testlet",
          function(.Object, id = NULL, item_list = NULL, model = "BTM",
                   parameters = NULL, se_parameters = NULL, content = NULL,
                   misc = NULL, ...) {
  .Object <- callNextMethod(.Object, ...)
  .Object@id <- id
  if (is(item_list, "Itempool")) .Object@item_list <- item_list else
    stop("'item_list' should be an 'Itempool' object." )
  if (is.character(model) && model %in% names(tmodels))
    .Object@model <- model else
      stop(paste0("Invalid model value. Model name should be one of the ",
                  'following:\n"', paste0(names(tmodels), collapse = '", "'),
                  '"'))
  .Object@model <- model
  .Object@parameters <- parameters
  .Object@se_parameters <- se_parameters
  .Object@content <- content
  .Object@misc <- misc

  # Check validity of the object
  validObject(.Object)
  .Object
})




############################################################################@###
############################# setValidity(Testlet) #############################
############################################################################@###
setValidity(
  Class = "Testlet",
  function(object)
  {
    # # ----------------------- Check item_list -------------------------------- #
    # # id's of all elements should be unique and they cannot be NULL.
    # # Itempool cannot be empty.
    # if (length(object@item_list) == 0)
    #   stop("Invalid 'item_list'. 'item_list' of Testlet object cannot be empty. ")
    # # All of the elements of the list should be "Item" class.
    # if (!all(sapply(object@item_list, FUN = "is.Item") ))
    #   stop(paste0("Invalid elements in item_list. All of the elements of
    #               'item_list' should be an 'Item' class object."))
    # # Item id's in item_list should be unique.
    # item_ids <- unlist(lapply(object@item_list, FUN = function(x) x@id))
    # if (is.null(item_ids) ||
    #       (length(unique(item_ids)) != length(object@item_list)))
    #   stop(paste0("Invalid Item id's. Each 'Item' object in the 'Testlet' ",
    #               "should have a unique id."))

    # ----------------------- Check id --------------------------------------- #
    # The length of testlet id should be 1 or NULL, and it cannot be NA
    if (!is.null(object@id) && (
      length(object@id) > 1 || is.na(object@id)))
      stop(paste0("Invalid testlet 'id'."))


    # ----------------------- Check model ------------------------------------ #
    # Check the model, currently only irt1PM, irt2PM,
    # irt3PM, irt4PM, mirt1PM, mirt2PM and mirt3PM is used
    if (is.null(object@model) || (length(object@model) != 1) ||
          !(object@model %in% names(tmodels)))
      stop(paste0("Invalid model name. Testlet model should be specified ",
                  "correctly. It can be either: ",
                  paste0(names(tmodels), collapse = ", ")))

    # ----------------------- Check parameters ------------------------------- #
    # Object parameters cannot be NULL or NA
    if (object@model != "BTM" &
        (is.null(object@parameters) || any(is.na(object@parameters))))
      stop(paste0("Invalid testlet parameter. Testlet parameters cannot be ",
                  "NULL or NA for models except 'BTM'."))

    # Number of parameters should be as specified in tmodels parameters
    if (length(object@parameters) != length(tmodels[[object@model]]$parameters))
      stop(paste0("Invalid testlet parameters. Number of parameters for ",
                  object@model," model should be ",
                  length(tmodels[[object@model]]$parameters), "."))

    # Check for proper naming of parameters. Parameter names should be unique
    # and all should correspond one-to-one with tmodels[[model_name]]$parameters
    if (object@model != "BTM") {
      if (is.null(parNames <- names(object@parameters)))
        stop(paste0("Invalid parameter names. Parameter names of Testlet class ",
                    "cannot be NULL except for 'BTM' model. Please give relevant names."))
      if ((length(parNames) !=  length(tmodels[[object@model]]$parameters)) ||
                   length(unique(parNames)) !=  length(tmodels[[object@model]]$parameters))
        stop(paste0("Invalid parameter names. Parameter names of Testlet class ",
                    "should be unique and complete. Please give relevant names."))
      if (!all(parNames %in% tmodels[[object@model]]$parameters))
        stop(paste0("Invalid parameter names. Parameter names for ",
                    object@model," model should be ",
                    paste0(tmodels[[object@model]]$parameters, collapse = ", "),
                    ". Please give relevant names."))
    }
  })
