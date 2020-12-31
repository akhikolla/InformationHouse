
###############################################################################@
############################# is.Testlet #######################################
###############################################################################@
#' Check whether an object is a \code{\link{Testlet-class}} object
#'
#' @param x an object that is checked for being a member of 'Testlet' class
#'
#' @export
#'
#' @rdname is
#'
#' @author Emre Gonulates
#'
is.Testlet <- function(x){is(x, "Testlet")}


###############################################################################@
############################# length (Testlet) #################################
###############################################################################@
#' Find the length of a \code{\link{Testlet-class}} object
#'
#' @export
#'
#' @rdname length
#'
setMethod(f = "length", signature = "Testlet",
          definition = function(x) length(x@item_list))


###############################################################################@
############################# testlet ##########################################
###############################################################################@
#' Creates a \code{\link{Testlet-class}} object
#'
#' @description Create a \code{\link{Testlet-class}} object. It is recommended
#'   to use this function to create new \code{\link{Testlet-class}} objects.
#'
#' @param ... The object that is desired to be converted to a \code{Testlet}
#'          object. Also additional arguments related to the \code{Testlet}.
#'
#' @return An \code{\link{Testlet-class}} object.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#' @importFrom methods validObject new
#'
#' @author Emre Gonulates
#'
#' @examples
#' ip <- itempool(a = c(1, 1.4), b = c(-2, 1))
#' testlet(ip, id = "T1")
#' testlet(ip, id = "T1", content = "Algebra")
#' # Add misc field to the testlet:
#' testlet(ip, id = "T1", misc = list(form = "A1", operational = TRUE,
#'                                    admin_date = as.Date("2020-08-01")))
#'
#' # Add misc field to the testlet items:
#' testlet(itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
#'                  misc = list(list(sympson_hetter_k = .8, form = "B1"),
#'                              list(sympson_hetter_k = .9))),
#'         id = "t1")
"testlet" <- function(...)
{
  # This function will coerce objects to "Testlet" class
  args <- list(...)
  x = args[[1]]
  if (is(x, "Itempool")) {
    # Assign default values for slots
    if (!is.null(names(args)) && "model" %in% names(args)) {
      if ((is.character(args$model)) && (length(args$model) == 1) &&
          (args$model %in% names(tmodels))) {
        model <- args$model
      } else if (!is.null(args$model)) {
        stop(paste0("Invalid Model specification. Model name should be one of ",
                    "following: ", paste0(names(tmodels), collapse = ", "), "."))
      } else model <- "BTM"
    } else model <- "BTM"
    # Add parameters, if exists
    if (!is.null(names(args)) && "parameters" %in% names(args))
      parameters <- args$parameters else parameters <- NULL
    # Add id, if exists
    if (!is.null(names(args)) && "id" %in% names(args)) id <- args$id else
       id <- NULL
    # Add content, if exists
    if (!is.null(names(args)) && "content" %in% names(args))
      content <- args$content else content <- NULL
    # Add misc, if exists
    if (!is.null(names(args)) && "misc" %in% names(args))
      misc <- args$misc else misc <- NULL
    # Add se_parameters, if exists
    if (!is.null(names(args)) && "se_parameters" %in% names(args))
      se_parameters <- args$se_parameters else se_parameters <- NULL
    return(new(Class = "Testlet", item_list = x, id = id, model = model,
               parameters = parameters, se_parameters = se_parameters,
               misc = misc, content = content))
  } else if (is(x, "Item")) {
    # Extract the elements that are item, create a list of them and assign as
    # the first element of arguments
    args[[1]] <- args[sapply(args, is.Item)]
    args[sapply(args, is.Item)] <- NULL
    return(do.call(testlet, args))
  } else if (is(x, "list")) {
    # All elements of the list should be 'Item' object
    if (all(sapply(x, is.Item))) {
      args[[1]] <- itempool(x)
      return(do.call(testlet, args))
    }
  } else if (inherits(x, c('matrix', 'data.frame'))) {
    args[[1]] <- itempool(x)
    return(do.call(testlet, args))
  }
  stop("Cannot coerce this object to a 'Testlet' object.")
}


###############################################################################@
############################# $ method #########################################
###############################################################################@
#' Access slots of a \code{\link{Testlet-class}} object
#'
#' @param x A \code{\link{Testlet-class}} object from which to extract
#'   element(s) or in which to replace element(s).
#'
#' @param name Name of the parameter.
#'          Available values:
#'          \describe{
#'            \item{\strong{\code{'id'}}}{Get the \code{id} of the testlet}
#'            \item{\strong{\code{'content'}}}{Get the \code{content} of
#'              the testlet.}
#'            \item{\strong{\code{'model'}}}{Get the \code{model} of the
#'              testlet.}
#'            \item{\strong{\code{'item_models'}}}{Get the \code{model} of the
#'              items within the testlet.}
#'            \item{\strong{\code{'parameters'}}}{Get the \code{parameters} of
#'              the testlet.}
#'            \item{\strong{\code{'se_parameters'}}}{Get the
#'              \code{se_parameters} of the testlets.}
#'            \item{\strong{\code{'item_list'}}}{Get the list of
#'              \code{\link{Item-class}} objects of the testlet. Returns a
#'              \code{list} object.}
#'            \item{\strong{\code{'max_score'}}}{Returns the maximum score
#'              obtainable by all of the items within the testlet.}
#'          }
#'
#' @return This operation will return the desired slot.
#'
#' @export
#' @examples
#' t1 <- testlet(generate_ip(n = 3), id = "my-testlet", content = "Algebra")
#' t1$model
#' t1$id
#' t1$item_list
#' t1$content
#' t1$item_models
setMethod("$", "Testlet",
          function(x, name)
          {
            switch(name,
                   'id'= return(x@id),
                   'item_list'= return(x@item_list@item_list),
                   'content' = return(x@content),
                   'model'= return(x@model),
                   'item_models'= return(x@item_list$model),
                   'parameters'= return(x@parameters),
                   'se_parameters'= return(x@se_parameters),
                   'max_score'= return(x@item_list$max_score),
                   x@parameters[[name]]
                   )
          })

###############################################################################@
############################# Subsetting 'Testlet' objects with "[" ############
###############################################################################@
#' Subset \code{\link{Testlet-class}} object
#'
#' @param x \code{Testlet} object from which to extract element(s) or in
#'   which to replace element(s).
#' @param i indices specifying elements to extract or replace.
#' @param j This will not be used in \code{\link{Testlet-class}} objects.
#' @param ... Parameters to be passed to the function.
#' @param drop (From R manual:) For matrices and arrays. If TRUE the result is
#' coerced to the lowest possible dimension (see the examples). This only works
#' for extracting elements, not for the replacement. See drop for further
#' details.
#'
#' @return An list object with elements from 'Item' class.
#'
#' @export
#'
#' @importFrom methods new
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
setMethod("[", c(x = "Testlet", i = "ANY", j = "missing", drop = "ANY"),
          function(x, i, j, ..., drop = TRUE)
          {
            return(new(Class = "Testlet", item_list = x@item_list[i],
                       id = x@id, model = x@model, parameters = x@parameters,
                       se_parameters = x@se_parameters, content = x@content))
          })


###############################################################################@
############################# Subsetting 'Testlet' objects with "[[" ###########
###############################################################################@
#' Access the items of a \code{\link{Testlet-class}} object.
#'
#' @param x A \code{Testlet}  object from which to extract element(s) or in
#'   which to replace element(s).
#' @param i indices specifying elements to extract or replace.
#' @param j This will not be used in \code{\link{Testlet-class}} objects.
#' @param ... Additional parameters to be passed to the function.
#'
#' @return An  object with elements from 'Item' class.
#'
#'
#' @export
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
setMethod("[[", c("Testlet", "numeric", "missing"),
          function(x, i, j, ...)
          {
            return(x@item_list[[i]])
          })


###############################################################################@
############################# Setting 'Testlet' objects with "[[<-"  ###########
###############################################################################@
#' This function sets the elements of a Testlet objects.
#'
#' @param x A \code{Testlet} object from which to extract element(s) or
#'   in which to replace element(s).
#' @param i indices specifying elements to extract or replace.
#' @param j This will not be used in \code{\link{Testlet-class}} objects.
#' @param value An \code{Item} object.
#' @param ... Additional parameters to be passed to the function.
#'
#' @return An \code{teslet} object with elements from 'Item' class.
#'
#'
#' @export
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
setMethod("[[<-", signature = c("Testlet", "numeric", "missing"),
          function(x, i, j, value)
          {
            x@item_list[[i]] <- value
            return(testlet(x@item_list, id = x@id, model = x@model,
                              parameters = x@parameters,
                              se_parameters = x@se_parameters,
                              content = x@content))
          })


###############################################################################@
############################# print.Testlet ####################################
###############################################################################@
#' Print a \code{\link{Testlet-class}} object
#'
#' @param x An \code{\link{Testlet-class}} object to be printed.
#' @param ... Passed parameters.
#'
#' @author Emre Gonulates
#'
#' @keywords internal
.print.Testlet <- function(x, ...)
{
  cat("An object of class 'Testlet'.\n")
  # if (is.null(x@id)) cat("Item id: - \n",sep = "") else
  if (!is.null(x@id))
    cat(paste0(format_text("ID:", italic = TRUE), "      ", x@id, "\n"))
    # cat("\033[3;mID:      \033[0;m", x@id, "\n",sep = "")

  # Color print: https://stackoverflow.com/a/57031762/2275286
  cat(paste0(format_text("Model:   ", italic = TRUE),
             format_text(x@model, bold = TRUE), "\n"))
  # cat("\033[3;mModel:  \033[0;m \033[1;m", x@model , "\033[0;m\n",sep = "")

  if (!is.null(x@content))
    cat(paste0(format_text("Content:  ", italic = TRUE),
               ifelse(length(x@content) == 1, paste0(x@content),
                      paste0(x@content, collapse = "; ")), "\n"))
    # cat("\033[3;mContent:\033[0;m ", ifelse(length(x@content) == 1,
    #                              paste0(x@content),
    #                              paste0(x@content, collapse = "; ")),
    #     "\n",sep = "")
  if (!is.null(x@parameters)) {
    cat(paste0(format_text("Model Parameters:", italic = TRUE), "\n"))
    # cat("\033[3;mModel Parameters:\033[0;m\n")
    pars <- unlist(x@parameters[sort(names(x@parameters))])
    print(pars)
    cat("\n")
    # cat("--------------------------")
    cat(paste0(rep('-',
                   # 26: the length of "An object of class 'Item'."
                   # 9: the length of maximum characters printed for a numeric obj
                   max(length(pars) * min(max(nchar(pars)), 9) + length(pars) - 1, 26)),
               collapse = ""))
    cat("\n")
  }
  if (!is.null(x@se_parameters)) {
    cat(paste0(format_text("Standard error of parameters:", italic = TRUE),
               "\n"))
    # cat("\033[3;mStandard error of parameters:\033[0;m \n",sep = "")
    # This needs to be fixed in case the order of names is not correct; and it
    # should show names.
    print(unlist(x@se_parameters))
  }
  if (!is.null(x@misc)) {
    cat(paste0(format_text("Misc:", italic = TRUE), "\n"))
    # cat("\033[3;mMisc:\033[0;m \n",sep = "")
    for (i in seq_along(length(x@misc))) {
      cat(paste0("  ", names(x@misc)[i], ": ", unlist(x@misc[i]), "\n"))
    }
  }
  cat(paste0("\n", format_text("Item List:", italic = TRUE), "\n"))
  # cat("\n\033[3;mItem List:\033[0;m\n")
  .print.Itempool(x@item_list, print_header = FALSE)
}


###############################################################################@
############################# show.Testlet #####################################
###############################################################################@
#' Show a \code{\link{Testlet-class}} object
#'
#' @param object A \code{\link{Testlet-class}} object that will be showed.
#'
#' @export
#'
#' @rdname show
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
setMethod("show", "Testlet", function(object) {.print.Testlet(object)})




