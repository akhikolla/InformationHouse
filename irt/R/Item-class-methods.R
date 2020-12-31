
###############################################################################@
############################# item #############################################
###############################################################################@
#' Create an \code{Item} object
#'
#' @description This function is used for creating \code{\link{Item-class}}
#'   objects.
#'
#' @param ... The item parameter arguments.
#' @param model The model that item \code{parameters} represents. Currently
#'   model can be:
#'   1PL, 2PL, 3PL, 4PL, M1PL, M2PL and M3PL, GRM, PCM or GPCM.
#'   Ideally, a model should be specified for the construction of an
#'   \code{\link{Item-class}} object.
#' @param id Item id. Default value is \code{NULL}.
#' @param parameters A list containing numeric vectors that represent item
#'   parameters. Depending on the model these can change.
#' @param se_parameters Standard error of item parameters.
#' @param content Content information for item.
#' @param misc This slot is a list where one can put any information about
#'  the item. For example, one can enter the id's of the enemies of the current
#'  item as \code{misc = list(enemies = c("i1", i2))}. Or, one can enter
#'  Sympson-Hetter exposure control parameter K:
#'  \code{misc = list(sympson_hetter_k = .75)}.
#'
#' @return An \code{Item} class object.
#'
#' @export
#'
#' @importFrom methods new
#'
#' @author Emre Gonulates
#'
#' @examples
#' # Create 2PL item:
#' item(a = 1.2, b = -0.94)
#' item(a = 1.2, b = -0.94, model = "2PL")
#' # Specify scaling constant D:
#' item(a = 1.2, b = -0.94, D = 1.7)
#'
#' # Add additional item specifications:
#' # Add id
#' item(a = 1.2, b = -0.94, id = "My-Item-1")
#' # Add content
#' item(a = 1.2, b = -0.94, id = "My-Item-1", content = "Geometry")
#' # Add additional parameter
#' item(a = 1.2, b = -0.94, misc = list(sympson_hetter_k = 1))
#' # Add any argument to 'misc' field
#' i1 <- item(a = 1.2, b = -0.94, id = "item1", content = "Earth Science",
#'            misc = list(key = "C", operational = TRUE, type = "MC",
#'                        enemies = c("i2", "i3")))
#' # Access fields
#' i1$misc
#' i1$misc$key
#' i1$misc$operational
#' i1$misc$enemies
#' i1$a
#' i1$b
#' i1$D
#' i1$parameters
#' i1$id
#' i1$content
#'
#'
#' # Rasch Model
#' item(b = 1.2)
#' item(b = 1.2, model = "Rasch")
#'
#' # 1PL model:
#' item(b = 1.2, model = "1PL")
#' item(b = 1.2, D = 1)
#'
#' # 3PL model:
#' item(a = 0.92, b = 2.7, c = 0.17)
#' item(a = 0.92, b = 2.7, c = 0.17, model = "3PL")
#' item(a = 0.92, b = 2.7, c = 0.17, D = 1.7, model = "3PL")
#'
#' # 4PL model:
#' item(a = 0.92, b = 2.7, c = 0.17, d = 0.98)
#' item(a = 0.92, b = 2.7, c = 0.17, d = 0.98, model = "4PL")
#' item(a = 0.92, b = 2.7, c = 0.17, d = 0.92, D = 1.7, model = "4PL")
#' item(parameters = list(a = 0.92, b = 2.7, c = 0.17, d = 0.92, D = 1.7),
#'      model = "4PL")
#'
#' # Create a GRM model
#' item(a = 1.9, b = c(-1, 0.82, 1.5), model = "GRM")
#' item(parameters = list(a = 1.9, b = c(-1, 2), D = 1), model = "GRM")
#'
#' # Create a GPCM model
#' item(a = 1.9, b = c(-1.6, -0.09, 1.25), model = "GPCM")
#' item(parameters = list(a = 1.9, b = c(-1, 2), D = 1), model = "GPCM")
#'
#' # Create a GPCM2 model (Reparametrized GPCM model)
#' item(a = 1.9, b = 0.65, d = c(-1.6, -0.09, 1.25), model = "GPCM2")
#' item(parameters = list(a = 1.9, b = 0.65, d = c(-1.6, -0.09, 1.25), D = 1.7),
#'      model = "GPCM2")
#'
#' # Create a PCM model
#' item(b = c(-0.7, 0.72, 1.9), model = "PCM")
#' item(parameters = list(b = c(-1, 2)), model = "PCM")
#'
#' # Add additional arguments to items
#' i1 <- item(a = 1.2, b = 2)
#' i1 <- item(i1, id = "new_item_id", content = "Algebra")
item <- function(..., model = NULL, id = NULL, parameters = NULL,
                 se_parameters = NULL, content = NULL, misc = NULL) {
  args <- list(...)
  clean_to_vector <- function(x) unname(unlist(as.vector(x)))

  # Check if the first element is already an Item object, then the other
  # slots can be updated.
  if (length(args) > 0 && is(args[[1]], "Item")) {
    object <-  args[[1]]
    # Check if the model is NULL.
    # If it is NULL, first check whether parameters argument exists, if it does,
    # try to use that to get the item model. If 'parameters' argument is NULL
    # also, try the '...' argument to get item parameters and the model.
  } else if (is.null(model)) {
    # Check parameters, if it exists use it to find the model name. If it is
    # also NULL then try to get the arguments '...' and match the model
    # that matches the most.
    if (is.null(parameters)) {
      if (is.null(names(args))) {
        stop(paste0("Insufficient parameters. Please give more information ",
                    "about the item parameters or specify the model ",
                    "explicitly."))
      }
      # First check whether there is b parameter, it means it is either UIRT
      # models or PCM, GRM or GPCM. Actually since GRM and GPCM has the same
      # item parameter names, if there is 'a' value and 'b' values then it
      # will be assumed that the item is "GRM". In order to create a GPCM,
      # user must specify 'model = "GPCM"'.
      if ("b" %in% names(args)) {
        # Check the length of "b" parameter, if it is 1, then it is a
        # dichotomous IRT model. Otherwise, it is a polytomous IRT model
        # Dichotomous item
        if (length(args$b) == 1 && is.numeric(clean_to_vector(args$b))) {
          # Check if "a" parameter exists
          if ("a" %in% names(args)) { # 2PL, 3PL or 4PL
            parameters <- list(a = clean_to_vector(args$a),
                               b = clean_to_vector(args$b))
            if ("c" %in% names(args)) { # 3PL or 4PL
              parameters$c <- clean_to_vector(args$c)
              if ("d" %in% names(args)) { # 4PL
                model <- "4PL"
                parameters$d <- clean_to_vector(args$d)
              } else { # 3PL
                model <- "3PL"
              }
            } else { # 2PL
              model <- "2PL"
            }
            # Check if D parameter exists:
            parameters$D <- ifelse("D" %in% names(args),
                                   clean_to_vector(args$D), default_D_scaling)
          } else if ("D" %in% names(args)) { # 1PL
            model <- "1PL"
            parameters <- list(b = clean_to_vector(args$b),
                               D = clean_to_vector(args$D))
          } else { # Rasch model
            model <- "Rasch"
            parameters <- list(b = clean_to_vector(args$b))
          }
        } else { # Polytomous item.
          # If there is "a" parameter then it is a "GRM" otherwise it is a
          # "PCM" item
          if ("a" %in% names(args)) { # A "GRM" item
            model <- "GRM"
            parameters <- list(
              a = clean_to_vector(args$a), b = clean_to_vector(args$b),
              D = ifelse("D" %in% names(args), clean_to_vector(args$D),
                         default_D_scaling))
          } else { # A "PCM" item
            model <- "PCM"
            parameters <- list(b = clean_to_vector(args$b))
          }
        }
      } else if ("d" %in% names(args)) { # It is MIRT model
        if ("a" %in% names(args)) { # M2PL or M3PL
          parameters <- list(a = clean_to_vector(args$a),
                             d = clean_to_vector(args$d))
          if ("c" %in% names(args)) { # M3PL
            model <- "M3PL"
            parameters$c <- clean_to_vector(args$c)
          } else { # M2PL
            model <- "M2PL"
          }
        } else { # M1PL
          model <- "M1PL"
          parameters <- list(d = clean_to_vector(args$d))
        }
        parameters$D <- ifelse("D" %in% names(args), clean_to_vector(args$D),
                               default_D_scaling)
      } else {
        stop(paste0("Insufficient parameters. Please give more information ",
                    "about the item parameters or specify the model ",
                    "explicitly."))
      }
      object <- new("Item", model = model, parameters = parameters)
    } else { # parameters is not NULL
      # Coerce parameters to "..." and rerun the function. The parameters in
      # "..." will be ignored. If there is multiple parameters in "..." and
      # "parameters" argument raise error.

      # Check whether "..." and "parameters" conflict
      if (is.null(names(parameters)) ||
          (!is.null(names(args)) && any(names(parameters) %in% names(args)))) {
        stop(paste0("Conflicting parameters in '...' and 'parameters' ",
                    "argument."))
      }
      # pass the function to the item function as argument
      args <- parameters
      return(do.call(item, c(args, model = model, id = id, parameters = NULL,
                             se_parameters = se_parameters, content = content,
                             misc = misc)))
    }

    # If there is a model and it is valid then create an item based on it.
  } else if (model %in% names(Pmodels)) {
    # Check whether 'parameters' argument is null or not. If it is NULL,
    # then look for the item parameters in the '...', except "D". If "D" is
    # missing by default, make it 1.
    # If it is not NULL, check whether it is valid. If 'parameters' argument
    # is valid, then use it in the model. If not raise an error.
    par_names <- names(Pmodels[[model]]$parameters)
    if (is.null(parameters)) {

      parameters <- list()

      par_sizes <- sapply(Pmodels[[model]]$parameters, function(i) i$size)

      error_message <- paste0(
          "Incomplete parameters. Make sure to provide parameters in the ",
          "following format: \nitem(model = '", model, "', ",
          paste0(par_names, collapse = " = , "), " = , .....)")

      for (i in 1:length(par_sizes)) {
        par_name <- names(par_sizes)[i]
        if (par_sizes[[i]] == 1) {
          if (is.null(names(args)) ||
              (!par_name %in% names(args) && par_name != "D")) {
            stop(error_message, call. = FALSE)
          } else if (par_name == "D" && par_name %in% names(args)) {
            parameters[[par_name]] <- clean_to_vector(args[[par_name]])
          } else
            parameters[[par_name]] <- clean_to_vector(args[[par_name]])
        } else if (par_sizes[[i]] > 1) {
          if (par_name %in% names(args)) {
            parameters[[par_name]] <- clean_to_vector(args[[par_name]])
          } else {
            jj <- NULL
            for (j in 1:par_sizes[[i]]) {
              if (all(paste0(par_name, 1:j) %in% names(args))) {
                jj <- j
              } else break
            }
            if (is.null(jj)) {
              stop(error_message, call. = FALSE)
            } else {
              parameters[[par_name]] <- unname(unlist(args[paste0(par_name, 1:jj)]))
            }
          }
        }
      }



      # if scaling parameter is necessary but it is not present, set a default
      # value to it
      if (is.null(parameters$D) && "D" %in% par_names)
        parameters$D <- default_D_scaling
      object <- new("Item", model = model, parameters = parameters)
    } else if (
      is.list(parameters) &&
      (
        # Full item parameters with D
        (length(parameters) == length(par_names) &&
         all(par_names %in% names(parameters))) ||
        # Full item parameters without D is fine, simply set D to default later
        ("D" %in% par_names &&
         length(parameters) == (length(par_names) - 1) &&
         all(setdiff(par_names, "D") %in% names(parameters)))
        )
      ) {
      # if scaling parameter is necessary but it is not present, set a default
      # value to it
      if (is.null(parameters$D) && "D" %in% par_names)
        parameters$D <- default_D_scaling
      object <- new("Item", model = model,
                    parameters = lapply(parameters, clean_to_vector))

      # Raise an error if "parameters" argument is not valid:
    } else {
       stop(paste0("Invalid item parameters. Item parameters for the model '",
                   model , "' should be a list with the following names: ",
                   paste0("'", par_names, "'", collapse = ", "), "."))
    }
    # If 'model' is not valid, raise an error.
  } else  {
    stop(paste0("Invalid 'model' specification. Model name should be one of ",
                "the following: ",
                paste0(names(Pmodels), collapse = ", "), "."))
  }

  # Add additional arguments
  # id
  if (!is.null(id)) {
    object@id <- id
  } else if ("id" %in% tolower(names(args)))
    for (id_name in c("ID", "Id", "iD"))
      if (id_name %in% names(args)) {
        object@id <- args[[id_name]]
        break
      }
  if (!is.null(content)) object@content <- content
  if (!is.null(misc)) object@misc <- misc
  if (!is.null(se_parameters)) object@se_parameters <- se_parameters
  if ("D" %in% names(args)) object@parameters$D <- args$D

  # Check the validity of the object and return the object if it is valid, else
  # raise an error.
  if (validObject(object))
    return(object)
}

###############################################################################@
############################# is.Item ##########################################
###############################################################################@
#' Check whether an object is an \code{\link{Item-class}}
#' @param x an object that is checked for whether being a
#'   \code{\link{Item-class}} object or not
#'
#' @export
#'
#' @importFrom methods is
#'
#' @rdname is
#'
#' @author Emre Gonulates
#'
#' @examples
#' i1 <- item(a = 1, b = 2)
#' is.Item(i1)
#' # Alternatively:
#' is(i1, "Item")
#'
#' # Not an item:
#' is.Item("abc")
#'
is.Item <- function(x){is(x,"Item")}


###############################################################################@
############################# print.Item #######################################
###############################################################################@
#' Print an \code{\link{Item-class}} object
#'
#' @param x An \code{\link{Item-class}} object to be printed.
#' @param ... Passed parameters.
#'
#' @author Emre Gonulates
#'
#' @keywords internal
#'
print.Item <- function(x, ...)
{
  cat("An object of class 'Item'.\n")
  # if (is.null(x@id)) cat("Item id: - \n",sep = "") else
  if (!is.null(x@id))
    cat(paste0(format_text("ID:", italic = TRUE), "      ", x@id, "\n"))
    # cat("\033[3;mID:      \033[0;m", x@id, "\n",sep = "")

  # Color print: https://stackoverflow.com/a/57031762/2275286
  cat(paste0(format_text("Model:", italic = TRUE), "   ",
             format_text(x@model, bold = TRUE), "\n"))
  # cat("\033[3;mModel:  \033[0;m \033[1;m", x@model , "\033[0;m\n",sep = "")

  if (!is.null(x@content))
    cat(paste0(format_text("Content:", italic = TRUE), " ",
               ifelse(length(x@content) == 1, paste0(x@content),
                      paste0(x@content, collapse = "; ")), "\n"))
    # cat("\033[3;mContent:\033[0;m ", ifelse(length(x@content) == 1,
    #                              paste0(x@content),
    #                              paste0(x@content, collapse = "; ")),
    #     "\n",sep = "")
  cat(paste0(format_text("Model Parameters:", italic = TRUE), "\n"))
  # cat("\033[3;mModel Parameters:\033[0;m\n")
  pars <- unlist(x@parameters[sort(names(x@parameters))])
  print(pars)
  cat("\n")

  if (!is.null(x@se_parameters)) {
    cat("\033[3;mStandard error of parameters:\033[0;m \n",sep = "")
    # This needs to be fixed in case the order of names is not correct; and it
    # should show names.
    print(unlist(x@se_parameters))
  }

  print_with_quotes <- function(x)
    ifelse(test = is(x, "character"), yes = paste0("\"", x, "\""),
           no = paste0(x))
  if (!is.null(x@misc)) {
    cat("\033[3;mMisc:\033[0;m \n",sep = "")
    for (i in seq_len(length(x@misc))) {
      if (length(x@misc[[i]]) == 1) {
        cat(paste0("  ", names(x@misc)[i], ": ",
                   print_with_quotes(unlist(x@misc[i])), "\n"))
      } else
        cat(paste0("  ", names(x@misc)[i], ": ",
                   paste0(sapply(unlist(x@misc[i]), print_with_quotes),
                          collapse = ", "), "\n"))
    }
  }
  # cat("--------------------------")
  cat(paste0(rep('-',
                 # 26: the length of "An object of class 'Item'."
                 # 9: the length of maximum characters printed for a numeric obj
                 max(length(pars) * min(max(nchar(pars)), 9) +
                       length(pars) - 1, 26)),
             collapse = ""), "\n")
}


###############################################################################@
############################# show.Item ########################################
###############################################################################@
#' Show an \code{\link{Item-class}} object
#'
#' @param object An \code{\link{Item-class}} object.
#'
#' @export
#'
#' @keywords internal
#'
#' @rdname show
#'
#' @author Emre Gonulates
#'
#' @importFrom methods show
#'
setMethod("show", "Item", function(object) {print.Item(object)})




###############################################################################@
############################# $ method #########################################
###############################################################################@
#' Get slots from an \code{\link{Item-class}} object.
#'
#' @param x An \code{\link{Item-class}} object.
#' @param name Name of the parameter.
#'
#'   Available values:
#'   \describe{
#'     \item{\strong{\code{'id'}}}{Extract \code{'id'} of an
#'       \code{\link{Item-class}} object.}
#'     \item{\strong{\code{'model'}}}{Extract the \code{'model'} of an
#'       \code{\link{Item-class}} object.}
#'     \item{\strong{\code{'parameters'}}}{Extract the \code{'parameters'} of an
#'       \code{\link{Item-class}} object.}
#'     \item{\strong{\code{'se_parameters'}}}{Extract the standard error of
#'       parameters of an \code{\link{Item-class}} object.}
#'     \item{\strong{\code{'content'}}}{Extract the \code{'content'} slot of an
#'       \code{\link{Item-class}} object.}
#'     \item{\strong{\code{'misc'}}}{Extract the \code{'misc'} slot of an
#'       \code{\link{Item-class}} object.}
#'     \item{\strong{\code{'max_score'}}}{Extract the maximum possible score
#'       of an \code{\link{Item-class}} object. Minimum score is assumed to be
#'       0.}
#'   }
#'
#' @return This operation will return the desired slot.
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @include Item-class.R
#'
#' @examples
#' item1 <- item(model =  "3PL", id = 'item23', content = 'Geometry',
#'               misc = list(enemies = c("item1", "item2")),
#'               parameters = list(b = 2, c = .12, a = 1.2, D = 1))
#' # Get individual parameters
#' item1$a
#' item1$b
#' item1$D
#' # Get item 'model'
#' item1$model
#' # Get all parameters
#' item1$parameters
#' # Get item 'id'
#' item1$id
#' # Get item content
#' item1$content
#' # Get misc values
#' item1$misc
#' # Get misc values
#' item1$max_score
setMethod("$", "Item",
          function(x, name)
          {
            switch(name,
                   "id" = return(x@id),
                   "content" = return(x@content),
                   "misc" = return(x@misc),
                   "model" = return(x@model),
                   "parameters" = return(x@parameters),
                   "se_parameters" = return(x@se_parameters),
                   "max_score" = return(get_max_possible_score_item_cpp(x)),
                   x@parameters[[name]]
                   )
          })


###############################################################################@
############################# $<- method #######################################
###############################################################################@
#' Set values to parameters or components of \code{\link{Item-class}} object
#'
#' @param x An \code{\link{Item-class}} object.
#' @param name Name of the parameter or component.
#' @param value The new value that will be assigned.
#'
#' @return This operation will not return anything.
#'
#' @export
#'
#' @importFrom methods new
#'
#' @author Emre Gonulates
#'
#' @include Item-class.R
#'
#' @examples
#' item <- new("Item", model =  "3PL", id = 'item23', content = 'Geometry',
#'             misc = list(enemies = c("item1", "item2")),
#'             parameters = list(b = 2, c = .12, a = 1.2, D = 1))
#' item$a <- 2
#' item$D <- 1.7
#' item$id <- "Itm-111"
#' item$content <- 'Algebra'
#' item$misc <- list(enemies = c("item5"))
setMethod(
  "$<-", "Item",
  function(x, name, value) {
    switch(name,
           'id'= {x@id <- value},
           'content' = {x@content <- value},
           'misc' = {x@misc <- value},
           # 'model'= x@model <- value,
           'parameters'= {
             par_names <- names(Pmodels[[x@model]]$parameters)
             if (is.list(value) &&
                 length(value) == length(par_names) &&
                 all(par_names %in% names(value))) {
               x@parameters <- value
               } else
                 warning(paste0(
                 "Invalid item parameters. Assignment cannot be made. ",
                 "Required item paramters are: ", paste0(
                   "'", par_names, "'", collapse = ", ")))
             },
           'se_parameters'= {x@se_parameters <- value},
           # The default is checking item
           {
             if (name %in% names(x@parameters))
               x@parameters[[name]] <- value
           }
           )
    if (validObject(x)) return(x)
  }
  )


