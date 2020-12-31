
###############################################################################@
############################# Itempool ########################################
###############################################################################@
#' Create an \code{Itempool} object
#'
#' @description This method creates a new \code{\link{Itempool-class}} object.
#'
#' @param ... The object that is desired to be converted to an 'Itempool'
#'          object. Also additional arguments related to the Itempool.
#'
#' @return An \code{\link{Itempool-class}} object.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#'
#' @importFrom methods validObject new
#'
#' @author Emre Gonulates
#'
#' @examples
#' # Create an item pool with two 2PL items
#' itempool(a = c(1, 1.4), b = c(-2, 1))
#' itempool(a = c(1, 1.4), b = c(-2, 1), model = "2PL")
#' # Set D parameter
#' itempool(a = c(1, 1.4), b = c(-2, 1), D = 1.7)
#' # Set item IDs
#' itempool(a = c(1, 1.4), b = c(-2, 1), id = c("i1", "i2"))
#' # Set content
#' itempool(a = c(1, 1.4), b = c(-2, 1), content = c("Algebra", "Geometry"))
#'
#' # Create GRM (Graded Response Model) items
#' # itempool(data.frame(a = rlnorm(10, 0, .3), b1 = rnorm(10), b2 = rnorm(10)),
#' #          model = "GRM")
#'
#' # Create a Rasch model item pool
#' itempool(b = c(-1, 0.2, 1.1), model = "Rasch")
#'
#' # Add 'misc' field:
#' ip <- itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
#'                 misc = list(list(sympson_hetter_k = .8),
#'                             list(sympson_hetter_k = .9)))
#' ip[[1]]  # First item of the item pool
#'
itempool <- function(...){
  # This function will create new "Itempool" objects.
  args <- list(...)
  # There are three possibilities to create an Itempool object:
  # (1) Give an argument like itempool(x = ....).
  # (2) Give an argument like itempool(<Object without Name>,....).
  # (3) Give arguments like item parameter names like
  #     itempool(a = runif(10, .5, 1.5), b = rnorm(2), ....).
  # If it looks like the third option, I'll create an 'x' object out of them.


  ## ------------------------  Functions ------------------------------------ ##
  # This function converts adds se_parameters given as a data frame or list
  # to a list of items (not testlets yet)
  #
  # @param item_list a list of items/testlets  (for now just items)
  # @param args The arguments vector supplied to the main function
  #
  # @return A list of items/testlets where the \code{se_parameters} slot
  #   updated based on the entry
  add_se_parameters <- function(item_list, args) {
    if (is.null(args$se_parameters)) return(item_list)
    se_parameters <- args$se_parameters
    # Get se_parameter name of each item and define acceptable names for
    # parameters
    se_par_names <- lapply(sapply(
      item_list, function(x) x@model),
      function(x) {
        p <- Pmodels[[x]]$parameters
        se <- names(p)[sapply(p, `[[`, "se")];
        return(c(se, paste0(se, "_se"), paste0("se_", se)))
        })
    # Convert matrix to data frame
    if (inherits(se_parameters, "matrix"))
      se_parameters <- as.data.frame(se_parameters)

    if (inherits(se_parameters, "data.frame")) {
      if (nrow(se_parameters) !=  length(item_list))
        stop("Invalid 'se_parameters. The number of rows of the ",
             "'se_parameters' should be ", length(item_list), ".")
      # Convert data frame rows to list elements
      se_parameters <- lapply(seq_len(nrow(se_parameters)),
           function(i) lapply(se_parameters, "[", i))
      # sanitize the names of the se parameters
      se_parameters <- lapply(se_parameters, function(x) {
        names(x) <- gsub("_se$|^se_", "", names(x));
        x
        })
    }
    # se_parameter should follow a certain structure
    if (!inherits(se_parameters, "list"))
      stop("Invalid 'se_parameters' provided to 'itempool()' function. ",
           "Please provide a list or data.frame for 'se_parameters'",
           call. = FALSE)
    item_list <- lapply(1:length(item_list), FUN = function(i) {
      item_list[[i]]@se_parameters <- se_parameters[[i]];
      return(item_list[[i]])})
    return(item_list)
  }

  ## -------------------------------------------------------------------------##
  # This function adds item id's given as a data frame or list
  # to a list of items (not testlets yet)
  #
  # @param item_list a list of items/testlets  (for now just items)
  # @param args The arguments vector supplied to the main function
  #
  # @return A list of items/testlets where the \code{id} slot
  #   updated based on the entry
  add_id <- function(item_list, args) {
    # ID can be supplied in two ways. The first, it can be supplied as an
    # argument to "args". Or it can be a column if the first element of the
    # function is a data.frame/matrix

    id <- NULL
    # Check whether 'id' exists in the arguments
    if ("id" %in% tolower(names(args))) {
      for (id_name in c("id", "ID", "Id", "iD"))
        if (id_name %in% names(args)) {
          id <- args[[id_name]]
          break
        }
    } else if (length(args) > 0 &&
               inherits(args[[1]], c('matrix', 'data.frame')) &&
               "id" %in% tolower(colnames(x))) {
      for (id_name in c("id", "ID", "Id", "iD"))
        if (id_name %in% colnames(args[[1]])) {
          id <- args[[1]][, id_name]
          break
        }
      if (is.null(id) && !is.null(rownames(args[[1]])))
        id <- rownames(args[[1]])
    }
    if (is.null(id)) {
      return(name_items(item_list))
    }
    id <- as.character(id)
    # Item id cannot be duplicated
    if (any(duplicated(id)))
      stop("Invalid item IDs. Item ID cannot be duplicated.", call. = FALSE)
    if (length(id) != length(item_list))
      stop(paste0("Invalid item IDs. Item ID length provided should be ",
                  length(item_list), "."), call. = FALSE)
    # Assign item ids
    item_list <- lapply(1:length(item_list), FUN = function(i) {
      item_list[[i]]@id <- id[i];
      return(item_list[[i]])})
    return(item_list)
  }

  ## -------------------------------------------------------------------------##
  # This function adds item content's given to a list of items (not testlets yet)
  #
  # @param item_list a list of items/testlets  (for now just items)
  # @param args The arguments vector supplied to the main function
  #
  # @return A list of items/testlets where the \code{content} slot
  #   updated based on the entry
  add_content <- function(item_list, args) {
    # content can be supplied in two ways. The first, it can be supplied as an
    # argument to "args". Or it can be a column if the first element of the
    # function is a data.frame/matrix

    content <- NULL
    # Check whether 'content' exists in the arguments
    if ("content" %in% tolower(names(args))) {
      content <- args[[names(args)["content" == tolower(names(args))][1]]]
    } else if (length(args) > 0 &&
               inherits(args[[1]], c('matrix', 'data.frame')) &&
               "content" %in% tolower(colnames(x))) {
      content <- args[[1]][, colnames(args[[1]])[
        "content" == tolower(colnames(args[[1]]))][1]]
    }
    if (is.null(content)) return(item_list)
    content <- as.character(content)
    # If the length of content is not equal to item_list, recycle content:
    content <- rep(content, length.out = length(item_list))

    # Assign item/testlet content
    item_list <- lapply(1:length(item_list), FUN = function(i) {
      item_list[[i]]@content <- content[i];
      return(item_list[[i]])})
    return(item_list)
  }

  ## -------------------------------------------------------------------------##
  # This function adds item 'misc' slot given to a list of items
  # (not testlets yet)
  #
  # @param item_list a list of items/testlets  (for now just items)
  # @param args The arguments vector supplied to the main function
  #
  # @return A list of items/testlets where the \code{misc} slot
  #   updated based on the entry
  add_misc <- function(item_list, args) {
    misc <- NULL
    # Check whether 'misc' exists in the arguments
    if ("misc" %in% tolower(names(args))) {
      misc <- args[[names(args)["misc" == tolower(names(args))][1]]]
    }
    if (is.null(misc)) return(item_list)

    # Make sure 'misc' is a list and all elements of the 'misc' are also list.
    if (!is.list(misc))
      stop("'misc' should be a list. If 'misc' values of items are ",
           "different than each other, then 'misc' should be a list of lists",
           "where each list element represents the 'misc' field of an item. ",
           "For example, if a 'time_limit' misc field is desired to be ",
           "assigned to each element of an item pool with two items, 'misc' ",
           "field should look like this: ",
           "'misc = list(list(time_limit = 120), list(time_limit = 180))'. ",
           "Whereas, if all items should have the same 'misc' field, the ",
           "following can be written: 'misc = list(test_form_id = \"F1\")'.",
           call. = TRUE)

    # If the length of misc is not equal to item_list, recycle content:
    if (!is.list(misc) || !all(sapply(misc, is.list)))
      misc <- list(misc)

    misc <- rep(misc, length.out = length(item_list))

    if (!is.list(misc) || !all(sapply(misc, is.list)) ||
        length(misc) !=  length(item_list))
      stop("Cannot add the misc field. Please provide a valid 'misc' value.")

    # Assign item/testlet content
    item_list <- lapply(1:length(item_list), FUN = function(i) {
      item_list[[i]]@misc <- misc[[i]];
      return(item_list[[i]])})
    return(item_list)
  }

  ## -------------------------------------------------------------------------##
  # This function adds D parameter to the items
  #
  # @param item_list a list of items/testlets  (for now just items)
  # @param args The arguments vector supplied to the main function
  #
  # @return A list of items/testlets where the appropriate "D" parameter added
  #   to the items
  add_D <- function(item_list, args) {
    D <- args$D
    if (is.null(D)) return(item_list)

    # If the length of D is not equal to item_list, recycle content:
    D <- rep(as.numeric(D), length.out = length(item_list))

    # Assign item/testlet content
    item_list <- lapply(1:length(item_list), FUN = function(i) {
      if ("D" %in% names(Pmodels[[item_list[[i]]@model]]$parameters))
        item_list[[i]]@parameters$D <- D[i];
      return(item_list[[i]])})
    return(item_list)
  }

  # Check whether the 'model' in function arguments is indeed a valid 'model'
  # name.
  check_model_argument <- function(model) {
    if (is.null(model) || !is.character(model) || length(model) != 1 ||
        (!model %in% names(Pmodels))) return(FALSE)
    return(TRUE)
  }



  ### --------------------- Start Main Function ---------------------------- ###

  # Scenario 1, first argument is an Itempool object and arguments are
  # relevant things to be added, like id, or D
  x <- args[[1]]
  if (inherits(x, 'Itempool')) {
    args[[1]] <- NULL
    # Add item ids
    x@item_list <- add_id(item_list = x@item_list, args = args)
    # Add content
    x@item_list <- add_content(item_list = x@item_list, args = args)
    # Add se_parameters
    x@item_list <- add_se_parameters(item_list = x@item_list, args = args)
    # Add misc
    x@item_list <- add_misc(item_list = x@item_list, args = args)
    # Add D
    x@item_list <- add_D(item_list = x@item_list, args = args)

    names(x@item_list) <- sapply(x@item_list, function(i) i@id)

    if (!all(sapply(x@item_list, FUN = "validObject")))
      stop("Invalid item. At least one 'Item' object is not valid." )
    if (!validObject(x))
      stop("Invalid 'Itempool'. At least one 'Item' object is not valid." )
    return(x)
  } else if (inherits(x, c("Item", "Testlet"))) {
    item_list_indicator <- sapply(args, function(m)
      is(m, "Item") | is(m, "Testlet"))
    item_list <- args[item_list_indicator]
    item_list <- name_items(item_list)
    if (length(item_list) == 1) {
      if (length(args) > 1)  x <- do.call("item", args)
      x <- name_items(list(x))
      names(x) <- sapply(x, function(i) i$id)
      # if (is.null(x@id))  x@id <- paste0("Item-", 1)
      result <- new("Itempool", item_list = x)
    } else if (length(item_list) > 1) {
      args[item_list_indicator] <- NULL
      args <- c(list(item_list), args)
      result <- do.call("itempool", args)
    }
    validObject(result)
    return(result)
  } else if (inherits(x, c("integer", "numeric"))) {
    # Check whether 'model' argument is in the arguments.
    # This part is only create items from Rasch, 1PL, 2PL, 3PL, 4PL models.
    #   -> (Yes): Check whether all of the item parameters for that model are
    #        in the argument:
    #        -> (Yes) Check:
    #             (1) whether the length of each element is the same and
    #             (2) whether each parameter set is numeric
    #             If there are parameters that are not
    #             Create itempool object.
    #        -> (No) Raise Error.
    #   -> (No): Check whether it is one of the Rasch, 1PL, 2PL, 3PL, 4PL models
    #        -> (Yes) Check
    #             (1) whether the length of each element is the same and
    #             (2) whether each parameter set is numeric
    #             Create itempool object.
    #        -> (No) Raise Error.

    # First check if model is in the arguments, if not, check whether the
    # arguments for Rasch, 1PL, 2PL, 3PL, 4PL models are complete. If they
    # are complete assign a model name, otherwise
    if (!"model" %in% names(args)) {
      for (m in c("4PL", "3PL", "2PL", "1PL")) {
        par_names <- names(Pmodels[[m]]$parameters)[
          sapply(Pmodels[[m]]$parameters, `[[`, "se")];
        if (all(par_names %in% names(args))) {
          args$model <- m
          break
        }
      }
    }
    if ("model" %in% names(args) && check_model_argument(args$model) &&
        args$model %in% names(Pmodels)[sapply(
          Pmodels, function(y) y$model_family == 'UIRT')]) {
      par_names <- names(Pmodels[[args$model]]$parameters)[
        sapply(Pmodels[[args$model]]$parameters, `[[`, "se")];
      # par_names <- Pmodels[[args$model]]$se_parameters
      if (all(par_names %in% names(args))) {
        # Check whether all of the parameter has the same length
        if (length(unique(sapply(args[par_names], length))) == 1 &&
            all(sapply(args[par_names], is.numeric))) {
          # Create a list of items using only the item parameters.
          n_items <- length(args[[par_names[1]]])
          result <- vector("list", n_items)
          pars <- do.call("cbind.data.frame", args[par_names])
          pars <- lapply(seq_len(nrow(pars)), function(i) lapply(pars, "[", i))
          result <- lapply(pars, function(i)
            do.call("item", args = c(model = args$model, i)))
          # Remove all parameter arguments from the arg
          args[par_names] <- NULL
          # Add results to the argument
          args <- c(list(name_items(result)), args)
          # Re-run the itempool() function with new arguments
          result <- do.call("itempool", args)
        } else {
          stop("All parameters supplied should be numeric and should have the ",
               "same length.", call. = FALSE)
        }
      } else {
        stop(paste0("Incomplete parameters. For '", args$model, "' model ",
                    "following parameters should be provided: ",
                    paste0("'", par_names, "'", collapse = ", "), "."))
      }
    } else
      stop(paste0("Itempool object cannot be created for '", args$model,
                  "' model. Please either choose a different model or ",
                  "use a different method to enter item parameter values.",
                  "see '?itempool'."), call. = FALSE)
    validObject(result)
    return(result)
  } else if (inherits(x, 'list')) {
    if (all(sapply(x, FUN = function(m) is(m, "Item") | is(m, "Testlet")))) {
      # Add id
      x <- add_id(item_list = x, args = args)
      # Add content
      x <- add_content(item_list = x, args = args)
      # Add se_parameters
      x <- add_se_parameters(item_list = x, args = args)
      # Add misc
      x <- add_misc(item_list = x, args = args)
      # Add D
      x <- add_D(item_list = x, args = args)

      names(x) <- sapply(x, function(i) i@id)
      result <- new("Itempool", item_list = x)
      validObject(result)
      return(result)
    } else stop("Invalid elements. All elements of the list should be an ",
                "'Item' or 'Testlet' object.")
  } else if (inherits(x, c('matrix', 'data.frame'))) {
    if (inherits(x, "matrix")) x <- as.data.frame(x)
    if ("model" %in% names(args)) {
      x <- cbind(x, model = args$model)
      args$model <- NULL
    }
    x <- lapply(seq_len(nrow(x)), function(i) lapply(x, "[", i))
    x <- lapply(x, function(i) do.call("item", i))
    args[[1]] <- x
    return(do.call("itempool", args))
  }
  stop("Cannot coerce this object to an 'Itempool' object.")
}

###############################################################################@
############################# name_items #######################################
###############################################################################@
#' Give Item class elements a unique id.
#' @description This function gives unique id's to elements of Item class. If
#'   there is no id's specified for an \code{\link{Item-class}} object, a
#'   default id will be given to that object. If some elements have id's
#'   already, uniqueness of the names will be checked. If they are not unique,
#'   an error will be issued. If some are unique and some are empty, a default
#'   id will be given to the empty ones.
#' @param item_list A list consist of \code{\link{Item-class}} or
#'   \code{\link{Testlet-class}} class objects.
#'
#' @return A list with \code{\link{Item-class}} object, which are all named.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @author Emre Gonulates
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' item1 <- item(a = 1.12, b = -2.1, c = 0.28)
#' item2 <- item(a = 2, b = 3.2, c = 0.21)
#' item4 <- item(a = 1.12, b = -1.23, c = .2, id = "I-21")
#' item5 <- item(a = 0.84, b = 2.23, c = .25, id = "I-22")
#'
#' item_list <- list(item1, item2)
#' name_items(item_list)
#'
#' item_list <- list(item1, item4, item2)
#' name_items(item_list)
#'
#' item_list <- list(item4, item5)
#' name_items(item_list)
#'
#' # Following code will give an error
#' # item_list <- list(item4, item4)
#' # name_items(item_list)
name_items <- function(item_list)
{
  # Stop if all elements of the list is not \code{\link{Item-class}} objects
  stopifnot(all(sapply(item_list,
                       FUN = function(x) is(x, "Item") | is(x, "Testlet"))))
  # Get item id's of only non-testlet items
  itemIDs <- unlist(sapply(
    item_list,
    FUN = function(itemObject) if (is(itemObject, "Item")) itemObject@id))
  # Check whether there is a testlet in the item_list:
  has_testlet <- any(sapply(item_list, is, "Testlet"))
  # If there is a testlet get id's of items within testlet.
  if (has_testlet)
    for (t in which(sapply(item_list, is, "Testlet")))
      itemIDs <- c(itemIDs, unlist(sapply(item_list[[t]]@item_list,
                                          FUN = function(itemObject) itemObject@id)))

  if (any(duplicated(itemIDs)))
    stop(paste0("\nInvalid id's. There are duplicated item id's. Correct ",
                "them before proceeding. Duplicated id's are: \n",
                paste0(itemIDs[duplicated(itemIDs)], collapse = ", ")),
         call. = FALSE)

  counter <- 1
  testlet_counter <- 1
  for (i in seq_len(length(item_list))) {
    if (is(item_list[[i]], "Item")) {
      while (is.null(item_list[[i]]@id)) {
        temp_id <- paste0("Item-", counter)
        if (!temp_id %in% itemIDs) {
          item_list[[i]]@id <- temp_id
          itemIDs <- c(itemIDs, temp_id)
        }
        counter <- counter + 1
      }
    } else if (is(item_list[[i]], "Testlet")) {
      # Give an id to testlet
      while (is.null(item_list[[i]]@id)) {
        temp_id <- paste0("Testlet-", testlet_counter)
        if (!temp_id %in% itemIDs) {
          item_list[[i]]@id <- temp_id
          itemIDs <- c(itemIDs, temp_id)
        }
        testlet_counter <- testlet_counter + 1
      }

      for (j in seq_len(length(item_list[[i]]))) {
        while (is.null(item_list[[i]]@item_list[[j]]@id)) {
          temp_id <- paste0("Item-", counter)
          if (!temp_id %in% itemIDs) {
            item_list[[i]]@item_list[[j]]@id <- temp_id
            itemIDs <- c(itemIDs, temp_id)
          }
          counter <- counter + 1
        }
      }
    }
  }

  # if (length(itemIDs) == 0)
  # {
  #   for (obj in item_list) {
  #     if (is(obj, "Item"))
  #   }
  #
  #
  #   itemIDs <- paste0("Item-", 1:length(item_list))
  #
  #   item_list <- sapply(seq_along(itemIDs),
  #                            FUN = function(a) {
  #                              item_list[[a]]@id <- itemIDs[a]
  #                              return(item_list[[a]])
  #                            })
  # } else if  {
  # } else if (length(itemIDs) < length(item_list))
  # {
  #   missingIDs <- which(sapply(sapply(
  #     item_list, FUN = function(itemObject) itemObject@id),
  #     FUN = "is.null"))
  #   itemIDs <- setdiff(paste0("Item-", 1:length(item_list)), itemIDs)
  #   item_list[missingIDs] <- sapply(
  #     seq_along(missingIDs), FUN = function(a) {
  #       item_list[[missingIDs[a]]]@id <- itemIDs[a]
  #       return(item_list[[missingIDs[a]]]) })
  # }
  return(item_list)
}



###############################################################################@
############################# concatenation of 'Item' objects ##################
###############################################################################@
.concatenate.Item.Itempool.obj <- function(x, ...) {
  args = list(x, ...)
  if (!all(sapply(args, inherits, c("Item", "Itempool", "Testlet"))))
    stop("All of the elements should be 'Item' class.")

  item_list = list()
  element_no = 0 # This designates the index of element at the final output ip
  for (i in 1:length(args)){
    if (is(args[[i]], "Item")) {
      element_no <- element_no + 1
      item_list[[element_no]] = args[[i]]
    } else if (is(args[[i]], "Itempool")) {
      for (j in 1:length(args[[i]])) {
        element_no <- element_no + 1
        item_list[[element_no]] <- args[[i]]@item_list[[j]]
      }
    } else if (is(args[[i]], "Testlet")) {
      element_no <- element_no + 1
      item_list[[element_no]] = args[[i]]
    }
  }
  # Name items in case they are missing.
  item_list <- name_items(item_list)
  names(item_list) <- sapply(item_list, function(i) i@id)
  return(methods::new(Class = "Itempool", item_list = item_list))
}

#' Concatenate \code{Item}, \code{Itempool} or \code{Testlet} objects and
#' return an Itempool object.
#'
#' If the elements do not have id fields, function will assign default names.
#'
#' @param x A list consist of \code{\link{Item-class}} objects.
#' @param ... Additional arguments
#'
#' @return An \code{\link{Itempool-class}} object.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @rdname c
#'
#' @method c Item
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @examples
#' item1 <- item(a = 1.12, b = -2.1, c = 0.28)
#' item2 <- item(a = 2, b = 3.2, c = 0.21)
#'
#' # Concatenate items
#' c(item1, item2)
#'
#' ip <- itempool(a = c(1, 1.2), b = c(1, 2), c = c(.2, .4))
#' # Concatenate items and an Itempool object
#' c(item1, ip)
#' c(item1, item2, ip)
#' c(ip, item1, item2)
setMethod("c", signature(x = "Item"), .concatenate.Item.Itempool.obj)

#' @rdname c
#' @method c Itempool
setMethod("c", signature(x = "Itempool"), .concatenate.Item.Itempool.obj)

#' @rdname c
#' @method c Testlet
setMethod("c", signature(x = "Testlet"), .concatenate.Item.Itempool.obj)

###############################################################################@
############################# Subsetting 'Itempool' objects with "[" ##########
###############################################################################@
#' Subset \code{Itempool} objects
#'
#' @param x An \code{\link{Itempool-class}} object from which to extract
#'   element(s) or in which to replace element(s).
#' @param i indices specifying elements to extract or replace.
#' @param j This will not be used in \code{\link{Itempool-class}} objects.
#' @param ... Parameters to be passed to the function.
#' @param drop (From R manual:) For matrices and arrays. If TRUE the result is
#' coerced to the lowest possible dimension (see the examples). This only works
#' for extracting elements, not for the replacement. See drop for further
#' details.
#'
#' @return An \code{\link{Itempool-class}} object with elements from
#'   \code{\link{Item-class}}.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
#' @examples
#' ip <- itempool(a = c(1.12, 2.1, 1.28), b = c(2, 3.2, 0.21),
#'                 id = c("i1", "i2", "i3"))
#'
#' ip[1]
#' # Create an Itempool using the first and third element:
#' ip[c(1, 3)] # Order is important
#' ip[c(3, 1)]
#' ip[-2]
#' ip[c(TRUE, FALSE, TRUE)]
#' ip[c("i2", "i1")]
#' # Recycle, i.e. get all elements
#' ip[TRUE]
setMethod("[", c(x = "Itempool", i = "ANY", j = "missing", drop = "ANY"),
          function(x, i, j, ..., drop = TRUE)
          {
            # If the length of i is the same as Itempool object and all of
            # its elements are character and none duplicated and all of them
            # are id's of elements (not Item id's within the testlets)
            if (is.character(i)) {
              if (!any(duplicated(i)) && all(i %in% x$id)) {
                x@item_list <- x@item_list[i]
              } else
                stop(paste0("Failed to subset using the given vector. Please ",
                            "use valid Item or testlet id's. There ",
                            "are either duplicated id's or non-existent id's ",
                            "in the subsetting vector provided."))
            } else if (is.logical(i) || is.numeric(i)) {
              # If there are NA values in the index then (1) if indices are
              # numeric, remove NAs. (2) if indices are logical, convert NAs to
              # FALSE
              if (any(is.na(i))) {
                if (is.numeric(i)) i <- i[!is.na(i)] else i[is.na(i)] <- FALSE
              }
              x@item_list <- x@item_list[i]
            }
            tryCatch({
              validObject(x)
              },
              error = function(e) {
                if (grepl("Item pool cannot be empty.", e$message))
                  stop("The selection did not match any Item/Testlet object.",
                       call. = FALSE)
                else stop(e$message, call. = FALSE)
              })
            return(x)
          })


###############################################################################@
############################# Subsetting 'Itempool' objects with "[[" #########
###############################################################################@
#' Subset \code{Itempool} objects
#'
#' @param x An \code{\link{Itempool-class}} object from which to extract
#'   element(s) or in which to replace element(s).
#' @param i indices specifying elements to extract or replace.
#' @param j This will not be used in \code{\link{Itempool-class}} objects.
#' @param ... Additional parameters to be passed to the function.
#'
#' @return An \code{\link{Item-class}} or \code{\link{Testlet-class}} object.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
#' @examples
#' item1 <- item(a = 1.12, b = -2.1, c = 0.28)
#' item2 <- item(a = 2, b = 3.2, c = 0.21)
#'
#' ip1 <- c(item1, item2)
#' ip1[[1]]
setMethod("[[", c("Itempool", "numeric", "missing"),
          function(x, i, j, ...)
          {
            result <- tryCatch({
              x@item_list[[i]]
              },
              error = function(e) {
                if (grepl("subscript out of bounds", e$message))
                  stop(paste0(
                    "Subscript out of bounds. Please use an index between ",
                    "1 and ", length(x), "."), call. = FALSE)
                NULL
              }
              )
            return(result)
          })


###############################################################################@
############################# Setting 'Itempool' objects with "[[<-"  #########
###############################################################################@
#' Set the elements of an \code{Itempool} objects.
#'
#' @param x \code{\link{Item-class}} object from which to extract element(s) or
#'   in which to replace element(s).
#' @param i indices specifying elements to extract or replace.
#' @param j This will not be used in \code{\link{Itempool-class}} objects.
#' @param value An \code{\link{Item-class}} object.
#' @param ... Additional parameters to be passed to the function.
#'
#' @return An \code{\link{Itempool-class}} object with elements from
#'   \code{\link{Item-class}} class.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
#' @examples
#' item1 <- item(a = 1.12, b = -2.1, c = 0.28)
#' item2 <- item(a = 2, b = 3.2, c = 0.21)
#' ip <- c(item1, item2)
#' item3 <- item(a = 1, b = -.2, c = 0.4)
#' ip[[2]] <- item3
setMethod("[[<-", signature = c("Itempool", "numeric", "missing"),
          function(x, i, j, value)
          {
            # Make sure the 'value' is either Item or Testlet
            if (!is(value, "Testlet") && !is(value, "Item"))
              stop(paste0("Invalid assignment. All elements of the list ",
                          "should be an 'Item' or 'Testlet' object."))
            x@item_list[[i]] <- value
            x <- itempool(x@item_list)
            validObject(x)
            return(x)
          })


###############################################################################@
############################# $<- method (Itempool) ###########################
###############################################################################@
#' Set values to parameters or components of 'Item' class.
#'
#' @param x \code{\link{Itempool-class}} object.
#' @param name Name of the parameter or component. Currently only \code{misc},
#'          \code{id}, \code{content}, \code{item_list} are available.
#' @param value The new value that will be assigned.
#'          \itemize{
#'            \item For \code{id}, the value should be a list of strings that
#'                  has the same length as the length of the
#'                  \code{\link{Itempool-class}} object. There should not be
#'                  any duplicated id's.
#'            \item For \code{content}, the value should be either \code{NULL}
#'                  or a list of strings that has the same length as the
#'                  length of the \code{\link{Itempool-class}} object.
#'            \item For \code{item_list}, the value should be a list of
#'                  \code{\link{Item-class}} or \code{\link{Testlet-class}}
#'                  objects.
#'            \item For \code{misc}, the value should be a list.
#'          }
#'
#'
#' @return This operation will return an \code{\link{Itempool-class}} object.
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
setMethod(
  "$<-", "Itempool",
  function(x, name, value) {
    switch(name,
           'id'= {
             # Check whether 'value' is character vector and the length of
             # the 'value' is equal to the length of the Itempool
             if (!is.character(value) || length(value) != length(x))
               stop(paste0("'id' should be a character vector with length ",
                           "equal to ", length(x), "."))
             # There shouldn't be a duplicated value.
             if (any(duplicated(value)))
               stop(paste0("'id' should not have any duplicated values. ",
                           "Check:\n", paste0("'", value[duplicated(value)],
                                              "'", collapse = ", ") ))
             for (i in seq_len(length(value)))
               x[[i]]@id <- value[i]
             },
           'content' = {
             # Check whether 'value' is either NULL or character vector and
             # the length of the 'value' is equal to the length of the Itempool
             if (!is.null(value) &&
                 (!is.character(value) || !length(value) %in% c(1, length(x))
                  ))
               stop(paste0("'content' should be a character vector with length",
                           " equal to ", length(x), "."))
             if (is.null(value)) {
               for (i in seq_len(length(x))) x[[i]]@content <- NULL
             } else if (length(value) == 1) {
               for (i in seq_len(length(x))) x[[i]]@content <- value
             } else
               for (i in seq_len(length(x))) x[[i]]@content <- value[i]
             },
           'item_list' = {
             if (!is(value, "Itempool") &&
                 all(sapply(value, function(i) is.Item(i) || is(i, "Testlet")))) {
               x@item_list <- value
             } else
               stop(paste0("'item_list' should be a list of 'Item' or ",
                           "'Testlet' objects."))
             },
           'misc' = {x@misc <- value},
           {
             model <- get_slot_itempool_cpp(ip = x, slotName = "model")
             if (all(model == model[1])) {
               result <- x$parameters
               if (name %in% colnames(result)) {
                 if (length(value) == 1) {
                   for (i in seq_len(length(x)))
                     x[[i]]@parameters[[name]] <- value
                 } else {
                   for (i in seq_len(length(x)))
                     x[[i]]@parameters[[name]] <- value[i]
                 }
               } else
                 stop(paste0("'", name, "' is not a valid name."))
             } else {
               stop(paste0("'", name, "' is not a valid name."))
             }
           }
           )
    validObject(x)
    return(x)
  }
  )



###############################################################################@
############################# $<- method (Testlet) #############################
###############################################################################@
#' Set values to parameters or components of 'Item' class.
#'
#' @param x A \code{\link{Testlet-class}} object.
#' @param name Name of the parameter or component.
#' @param value The new value that will be assigned.
#'
#' @return This operation will not return anything.
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
setMethod(
  "$<-", "Testlet",
  function(x, name, value) {
    switch(name,
           'id'= {x@id <- value; return(x)},
           'content' = {x@content <- value; return(x)},
           'misc' = {x@misc <- value; return(x)},
           'item_list' = {
             if (is(value, "Itempool")) {
               x@item_list <- value
               return(x)
             }
             },
           # 'model'= x@model <- value,
           'parameters'= {x@parameters <- value; return(x)},
           'se_parameters'= {x@se_parameters <- value; return(x)}, {
             if (name %in% names(x@parameters))
               x@parameters[[name]] <- value
             return(x)
           }
           )
  }
  )


###############################################################################@
############################# as.list (Itempool) ##############################
###############################################################################@
#' This function converts Itempool objects to a list object
#'
#' @param x an \code{\link{Itempool-class}} to be coerced to a list object
#' @param ... Additional parameters to be passed to the function.
#'
#' @return A list object with elements from 'Item' class.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @examples
#' item1 <- item(a = 1.12, b = -2.1, c = 0.28)
#' item2 <- item(a = 2, b = 3.2, c = 0.21)
#'
#' ip1 <- c(item1, item2)
#' as.list(ip1)
#'
as.list.Itempool <- function(x, ...) return(x@item_list)


###############################################################################@
############################# as.data.frame (Itempool) ########################
###############################################################################@
#' Convert an \code{\link{Itempool-class}} object into a \code{data.frame}.
#'
#' @description  This function converts \code{\link{Itempool-class}} objects to a
#'   \code{data.frame} object. All of the items in the item pool should be
#'   the same model. If the items belongs to different psychometric models,
#'   then, a \code{list} of items will be returned instead.
#'
#' @param x An \code{\link{Itempool-class}} object
#' @param row.names \code{NULL} or a character vector giving the row names for
#'   the data frame. Missing values are not allowed.
#' @param optional logical. If \code{TRUE}, setting row names and converting
#'   column names
#' @param ... additional arguments
#' @param include_se If \code{TRUE}, and items have \code{se_parameters},
#'   those will be included in the data frame.
#'
#' @return A data frame of items within each row. If all items cannot be
#'   coerced to a \code{data.frame}, an list of items will be returned and a
#'   warning will be raised.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @rdname as.data.frame
#'
#' @examples
#' ip1 <- generate_ip()
#' as.data.frame(ip1)
#'
#' ip2 <- generate_ip(n = 10, model = "GRM",
#'                    content = sample(c("G", "A"), 10, TRUE),
#'                    id = paste0("grm-i-", 1:10))
#' as.data.frame(ip2)
#'
#' t1 <- generate_testlet(n = 3, item_id_preamble = "t1")
#' t2 <- generate_testlet(n = 2, item_id_preamble = "t2")
#' ip3 <- c(ip1, t1, t2)
#' as.data.frame(ip3)
#'
#' ip4 <- c(ip2, ip3)
#' as.data.frame(ip4)
#'
#'
#' item1 <- item(a = 1.12, b = -2.1, c = 0.28)
#' item2 <- item(a = 2, b = 3.2, c = 0.21)
#'
#' ip1 <- c(item1, item2)
#' as.data.frame(ip1)
as.data.frame.Itempool <- function(x, row.names = NULL, optional = FALSE, ...,
                                   include_se = TRUE)
{
  args <- list(...)

  item_models <- x$resp_model

  # If Itempool object don't have any testlets and all item models are the same:
  if (length(unique(item_models)) == 1 &&
      !any(sapply(x@item_list, is, "Testlet"))) {
    pars <- x$parameters
    if (is.matrix(pars)) {
      result <- data.frame(pars, row.names = NULL, stringsAsFactors = FALSE)
      # This means that the models are the same but items are like GPCM or GRM
      # with different threshold number or 'mirt' with different number of
      # dimensions (i.e. a parameters)
    } else if (is.list(pars)) {
      result <- x$parameters
      # This step is to get the correct order of parameters.
      temp_df <- result[names(sort(sapply(result, length), decreasing = TRUE))]
      temp_df <- unique(unlist(sapply(temp_df, names)))
      temp_df <- data.frame(matrix(vector(), length(x), length(temp_df),
                                   dimnames = list(NULL, temp_df)))
      for (i in seq_len(length(x)))
        temp_df[i, names(result[[i]])] <- result[[i]]
      result <- temp_df
    }
    result <- cbind(id = x$id, model = unname(item_models), result)
  } else { # There are testlets in the item pool and there maybe multiple models
    temp_df <- data.frame(id = x$resp_id, testlet = NA, model = item_models,
                          stringsAsFactors = FALSE, row.names = NULL)
    dfs <- lapply(x@item_list, as.data.frame, include_se = FALSE)
    col_names <- setdiff(sort(unique(unlist(sapply(dfs, names)))),
                         c("id", "content", "testlet", "model"))
    result <- cbind(temp_df, matrix(vector(), nrow(temp_df), length(col_names),
                                     dimnames = list(NULL, col_names)))
    counter <- 1
    for (i in seq_len(length(dfs))) {
      result[counter:(counter + nrow(dfs[[i]]) - 1), names(dfs[[i]])] <- dfs[[i]]
      counter <- counter + nrow(dfs[[i]])
    }
    # Remove testlet column if there are no testlets:
    if (all(is.na(result$testlet))) result$testlet <- NULL
  }

  ### se_parameters ###
  se_pars <- x$se_parameters
  if (!is.null(se_pars) && include_se) {
    old_colnames <- colnames(result)
    colnames(se_pars) <- paste0(colnames(se_pars), "_se")
    result <- cbind(result, se_pars)
    # Create a new column order where the standard error of each parameter
    # follows the parameter:
    new_colnames <- c()
    for (i in old_colnames) {
      new_colnames <- c(new_colnames, i)
      candidate <- paste0(i, "_se")
      if (candidate %in% colnames(result)) {
        new_colnames <- c(new_colnames, candidate)
      }
    }
    result <- result[, new_colnames]
  }

  ### content ###
  content <- x$content
  if (!is.null(unlist(content)) &&
      (is.character(content) && (length(content) == nrow(result))) &&
      !"content" %in% colnames(result))
    result <- cbind.data.frame(result, content = unname(content),
                               stringsAsFactors = FALSE)

  ### misc ###
  # Check whether at least one item has 'misc' field
  misc <- x$resp_misc
  if (!is.null(misc)) {
    # Get all of the names of the unique fields
    for (misc_field in unique(unlist(sapply(misc, names)))) {
      temp <- lapply(misc, `[[`, misc_field)
      temp[sapply(temp, is.null)] <- NA
      # Only add columns which has length 1
      if (all(sapply(temp, length) == 1)) {
        result[[misc_field]] <- unlist(temp)
      }
    }
  }
  return(result)
}

###############################################################################@
############################# as.data.frame (Item) #############################
###############################################################################@
#' Convert an \code{\link{Item-class}} object into a \code{data.frame}.
#'
#' @description  This function converts \code{\link{Item-class}} objects to a
#'   \code{data.frame} object.
#'
#' @param x An \code{\link{Item-class}} object
#' @param row.names \code{NULL} or a character vector giving the row name for
#'   the data frame. Missing values are not allowed.
#' @param optional logical. If \code{TRUE}, setting row names and converting
#'   column names
#' @param ... additional arguments
#' @param include_se If \code{TRUE}, and items have \code{se_parameters},
#'   those will be included in the data frame.
#'
#' @return A data frame representation of the item.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @rdname as.data.frame
#'
#' @examples
#' item1 <- generate_item()
#' as.data.frame(item1)
as.data.frame.Item <- function(x, row.names = NULL, optional = FALSE, ...,
                               include_se = TRUE)
{
  return(as.data.frame(itempool(x), include_se = include_se))
}


###############################################################################@
############################# as.data.frame (Testlet) ##########################
###############################################################################@
#' Convert an \code{\link{Testlet-class}} object into a \code{data.frame}.
#'
#' @description  This function converts \code{\link{Testlet-class}} objects to a
#'   \code{data.frame} object. If testlet has an ID, an additional column
#'   will be created for the testlet ID.
#'
#' @param x An \code{\link{Testlet-class}} object
#' @param row.names \code{NULL} or a character vector giving the row name for
#'   the data frame. Missing values are not allowed.
#' @param optional logical. If \code{TRUE}, setting row names and converting
#'   column names
#' @param ... additional arguments
#' @param include_se If \code{TRUE}, and items have \code{se_parameters},
#'   those will be included in the data frame.
#'
#' @return A data frame representation of the item.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @rdname as.data.frame
#'
#' @examples
#' testlet1 <- generate_testlet()
#' as.data.frame(testlet1)
#' testlet2 <- generate_testlet(id = "T1")
#' as.data.frame(testlet2)
as.data.frame.Testlet <- function(x, row.names = NULL, optional = FALSE, ...,
                                  include_se = TRUE)
{
  output <- as.data.frame(x@item_list, include_se = include_se)
  if (!is.null(x@id))
    output <- cbind(testlet = x@id, output)
  return(output)
}

###############################################################################@
############################# is.ItemPool ######################################
###############################################################################@
#' Check whether an object is an  \code{\link{Itempool-class}} object
#'
#' @param x an object that is checked for being a member of 'Itempool' class
#'
#' @export
#'
#' @rdname is
#'
#' @author Emre Gonulates
#'
is.Itempool <- function(x){is(x,"Itempool")}


###############################################################################@
############################# length (Itempool) ###############################
###############################################################################@
#' Find the length of an \code{\link{Itempool-class}} object
#'
#' @param x an \code{\link{Itempool-class}} object
#'
#' @export
#'
#' @rdname length
#'
#' @author Emre Gonulates
#'
setMethod(f = "length", signature = "Itempool",
          definition = function(x) length(x@item_list))



###############################################################################@
############################# $ method #########################################
###############################################################################@
#' Get slots of the an \code{\link{Item-class}} object.
#'
#' @param x An \code{\link{Itempool-class}} object from which to extract
#'   element(s) or in which to replace element(s).
#' @param name Name of the parameter.
#'          Available values:
#'          \describe{
#'            \item{\strong{\code{'id'}}}{Extract \code{id}'s of all items
#'              and testlets.
#'              This will not extract the \code{id}'s of items within the
#'              testlet.}
#'            \item{\strong{\code{'content'}}}{Extract \code{content}'s of
#'              all items and testlets.
#'              This will not extract the \code{content}'s of items within the
#'              testlet.}
#'            \item{\strong{\code{'model'}}}{Extract \code{model}'s of
#'              all items and testlets.
#'              This will not extract the \code{model}'s of items within the
#'              testlet.}
#'            \item{\strong{\code{'misc'}}}{Extract \code{misc} parameters of
#'              all items and testlets.
#'              This will not extract the \code{misc} parameters of items
#'              within the testlet.}
#'            \item{\strong{\code{'item_list'}}}{Extract individual elements of
#'              item pool. If there are testlets in the item pool, a testlet
#'              will be an item of the resulting list. If individual items
#'              within the testlet is desired to be elements of the list, then
#'              use \code{$items}.}
#'            \item{\strong{\code{'items'}}}{Extract individual items
#'              within the item pool. If there are testlets in the item pool
#'              individual elements of the testlet will be extracted. Resulting
#'              list will only consist of \code{\link{Item-class}} objects.
#'            }
#'            \item{\strong{\code{'parameters'}}}{Extract \code{parameters}'s of
#'              all items and testlets.
#'              This will not extract the \code{parameters}'s of items within
#'              the testlet.}
#'            \item{\strong{\code{'se_parameters'}}}{Extract
#'              \code{se_parameters}'s of all items and testlets.
#'              This will not extract the \code{se_parameters}'s of items
#'              within the testlet.}
#'            \item{\strong{\code{'n'}}}{Return a list with three objects:
#'              \code{elements} the number of standalone items and testlets.
#'              \code{testlets} the number of Testlet objects.
#'              \code{items} the sum of the number of items within testlets and
#'              standalone items.
#'              }
#'            \item{\strong{\code{'max_score'}}}{Returns the maximum possible
#'              raw score of the item pool.
#'              }
#'            \item{\strong{\code{'resp_id'}}}{Extract
#'              \code{id}'s of all standalone items and items within the
#'              testlets. It will not return testlet \code{id}'s. This is
#'              mainly to get the \code{id}'s of items which has a response.
#'              }
#'            \item{\strong{\code{'resp_content'}}}{Extract
#'              \code{content}'s of all standalone items and items within the
#'              testlets. It will not return testlet \code{content}'s. This
#'              is mainly to get the \code{content}'s of items which has a
#'              response.
#'              }
#'            \item{\strong{\code{'resp_model'}}}{Extract
#'              \code{model}'s of all standalone items and items within the
#'              testlets. It will not return testlet \code{model}'s. This is
#'              mainly to get the \code{model}'s of items which has a response.
#'              }
#'            \item{\strong{\code{'resp_misc'}}}{Extract
#'              \code{misc} fields of all standalone items and items within
#'              the testlets. It will not return testlet \code{misc} fields.
#'              }
#'            \item{\strong{\code{'resp_item_list'}}}{Combine items that are
#'              not in a testlet and items within a testlet and return a list
#'              object. This list does not contain any \code{Testlet} objects.
#'              All of the elements are \code{Item} objects. If there are no
#'              testlets in the item pool, then this argument will be the
#'              same as \code{$item_list}.
#'              }
#'            \item{\strong{\code{'resp_max_score'}}}{Extract the maximum score
#'              each standalone item can get.
#'              }
#'          }
#'
#' @return This operation will return a numeric object.
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @include Item-class.R
#'
#' @examples
#' item1 <- methods::new("Item", model =  "3PL", id = 'item23',
#'                       content = 'Geometry',
#'                       parameters = list(b = 2, c = .12, a = 1.2, D = 1))
#' item1$a
#' item1$D
#' item1$model
#' item1$id
#' item1$content
setMethod(
  "$", "Itempool",
  function(x, name) {
    switch(
      name,
      'id'= return(get_slot_itempool_cpp(ip = x, slotName = "id")),
      'content' = {
        result <- get_slot_itempool_cpp(ip = x, slotName = "content")
        if (!is.null(result))
          result <- stats::setNames(
            result, get_slot_itempool_cpp(ip = x, slotName = "id"))
        return(result)
        },
      'model'= {
        return(stats::setNames(get_slot_itempool_cpp(ip = x, slotName = "model"),
                        get_slot_itempool_cpp(ip = x, slotName = "id")))
        },
      'misc'= {
        result <- lapply(x@item_list, function(i) i@misc)
        if (all(sapply(result, is.null))) result <- NULL
        return(result)
        },
      'item_list'= return(x@item_list),
      'items' = return(flatten_itempool_cpp(x)),
      'parameters'= {
        models <- x$model
        if (all(models == models[1])) {
          if (models[1] %in% names(Pmodels)[sapply(
            Pmodels, function(y) y$model_family == 'UIRT')]) {
            return(get_parameters_itempool_cpp(ip = x))
          } else if (models[1] %in% c(names(Pmodels)[sapply(
            Pmodels, function(y) y$model_family %in% c("PIRT", "MIRT"))])) {
            result <- sapply(1:length(x@item_list), FUN = function(i)
              unlist(x[[i]]@parameters))
            if (is.matrix(result)) {
              result <- t(result)
            } else if (is.atomic(result)) {
              result <- as.matrix(result)
              # par_names <- names(x@item_list[[1]]@parameters)
              colnames(result) <- names(x@item_list[[1]]@parameters)
              # This happens when for example items are polytomous items with
              # different number or thresholds or MIRT items with different
              # number of dimensions.
            } else if (is.list(result)) {
              names(result) <- x$id
              return(result)
            } else {
              result <- t(result)
            }
            rownames(result) <- get_slot_itempool_cpp(ip = x, slotName = "id")
          } else stop("This model is not implemented yet.")
        } else
        {
          result <- lapply(x@item_list, FUN = function(y) y@parameters)
          names(result) <- get_slot_itempool_cpp(ip = x, slotName = "id")
        }
        return(result)
        },
      'se_parameters'= {
        # Check whether there are any NULL SEs. If there are, temporarily,
        # set SE values for those items to NA
        null_se <- sapply(x, function(i) is.null(i@se_parameters))
        if (all(null_se)) {
          return(NULL)
        } else if (any(null_se)) {
          for (i in which(sapply(x, function(i) is.null(i@se_parameters)))) {
            # se_pars <- Pmodels[[x[[i]]@model]]$se_parameters
            se_pars <- names(Pmodels[[x[[i]]@model]]$parameters)[
              sapply(Pmodels[[x[[i]]@model]]$parameters, `[[`, "se")];


            x[[i]]@se_parameters <- stats::setNames(
			  rep(list(NA), length(se_pars)), se_pars)
          }
        }
        result <- sapply(x, function(i) unlist(i@se_parameters))
        # result <- sapply(x, function(i) as.data.frame(t(unlist(i@se_parameters))))
        if (inherits(result, "matrix")) {
           return(as.data.frame(t(result)))
        } else {
          result <- sapply(result, function(i) as.data.frame(t(i)))
          return(Reduce(rbind, Map(function(x) {
            x[, setdiff(unique(unlist(lapply(result, colnames))),
                        names(x))] <- NA;
            return(x)}, result)))
        }
        },
      'max_score'= return(sum(x$resp_max_score)),
      'resp_id'= {
        ids <- as.list(get_slot_itempool_cpp(ip = x, slotName = "id"))
        for (j in which(sapply(x, is, "Testlet"))) ids[[j]] <- x[[j]]@item_list$id
        return(unlist(ids))
        },
      'n'= {
        ip_size <- get_itempool_size(x)
        return(list(elements = unname(ip_size["elements"]),
                    testlets = unname(ip_size["testlets"]),
                    items = unname(ip_size["items"])))
        },
      'resp_content'= {
        content <- as.list(get_slot_itempool_cpp(ip = x, slotName = "content"))
        # Add content of the testlets. Here, if items of the testlet do not
        # have any content assigned, then they will be represented as NAs.
        for (j in which(sapply(x, is, "Testlet"))) {
          temp_content <- x[[j]]@item_list$content
          if (is.null(temp_content)) {
            content[[j]] <- rep(NA, length(x[[j]]@item_list))
          } else content[[j]] <- temp_content
        }
        return(stats::setNames(unlist(content), x$resp_id))
        },
      'resp_model'= {
        models <- as.list(get_slot_itempool_cpp(ip = x, slotName = "model"))
        for (j in which(sapply(x, is, "Testlet")))
          models[[j]] <- x[[j]]@item_list$model
        return(stats::setNames(unlist(models), x$resp_id))
        },
      'resp_misc'= {
        misc <- list()
        for (i in seq_along(x@item_list)) {
          if (is(x@item_list[[i]], "Testlet")) {
            temp_misc <- x@item_list[[i]]@item_list$misc
            if (is.null(temp_misc))
              temp_misc <- stats::setNames(
                vector("list", length(x@item_list[[i]]@item_list)),
                x@item_list[[i]]@item_list$id)
          } else {
            temp_misc <- stats::setNames(list(x@item_list[[i]]$misc),
                                         x@item_list[[i]]@id)
          }
          misc <- c(misc, temp_misc)
        }
        if (all(sapply(misc, is.null))) return(NULL)
        return(misc)
        },
      'resp_item_list'= {
        return(flatten_itempool_cpp(x))
        },
      'resp_max_score'= {
        return(get_max_possible_score_itempool_cpp(x))
        },
      # The default is checking whether individual parameters and extracting
      # them if all of the models are equal.
      {
        result <- as.data.frame(x)
        if (name %in% colnames(result)) {
          return(stats::setNames(result[, name], x$resp_id))
        } else return(NULL)
      }
    )
  }
  )


###############################################################################@
############################# .print.Itempool ##################################
###############################################################################@
# Prints an \code{\link{Itempool-class}} object
#
# @keywords internal
#
# @param x An \code{\link{Itempool-class}} object to be printed.
# @param ... Additional arguments.
# @param n maximum number of items to print. Default is \code{NULL}, where
#   all items are printed if the number of items are smaller than 20, otherwise
#   only first 10 items are printed.
# @param print_header Whether to print the object class in the first line.
#
# @author Emre Gonulates
# @noRd
# @examples
# (ip <- generate_ip(model = "3PL", n = 20))
# (ip <- c(generate_ip(model = "3PL", n = 20),
#          generate_testlet(id = "t1"), generate_testlet(id = "t2")))
.print.Itempool <- function(x, ..., n = NULL, print_header = TRUE) {
  # Make sure the printed object is valid
  tryCatch(validObject(x),
           error = function(e)
             stop(paste0("This is not a valid 'ItemPool' object:\n", e$message),
                  call. = FALSE))
  if (print_header)
    cat("An object of class 'Itempool'.\n")

  result <- as.data.frame(x)

  # Check ID and move it to row name
  # rownames(result) <- result$id
  # result$id <- NULL
  # Check whether all item models are the same
  if (length(unique(result$model)) == 1) {
    cat("Model of items: ", format_text(result$model[1], bold = TRUE), "\n",
        sep = "")
    # cat("Model of items: \033[1;m", result$model[1] , "\033[0;m\n",sep = "")
    result$model <- NULL
  }
  # Check if all D parameters are the same, if yes remove them.
  if (!is.null(result$D) && length(unique(result$D)) == 1) {
    cat(paste0("D = ", result$D[1],"\n"))
    result$D <- NULL
  }

  # Check if all content are the same and, if yes remove them.
  if (!is.null(result$content) && length(unique(result$content)) == 1) {
    cat(paste0("Content = ", result$content[1],"\n"))
    result$content <- NULL
  }

  # Print common fields
  if (nrow(result) > 1) {
    remove_cols <- c()
    for (i in seq_len(ncol(result)))
      if (length(unique(result[, i])) == 1) {
        cat(paste0(colnames(result)[i], " = ", result[1, i], "\n"))
        remove_cols <- c(remove_cols, i)
      }
    result[, remove_cols] <- NULL
  }

  cat("\n")
  has_testlet <- any(sapply(x@item_list, is, "Testlet"))
  print_ip_size <- function(x) {
    if (has_testlet) {
      n_testlets <- x$n$testlets
      n_standalone_items <- x$n$elements - x$n$testlets
      n_all_items <- x$n$items
      n_digits <- nchar(as.character(n_all_items))

      text_after_df <- sprintf(paste0("%32s = %", n_digits,"d\n"),
                               "Number of Testlets", n_testlets)
      text_after_df <- paste0(
        text_after_df,
        sprintf(paste0("%32s = %", n_digits,"d\n"),
                "Number of items within Testlets",
                n_all_items - n_standalone_items))
      text_after_df <- paste0(
        text_after_df,
        sprintf(paste0("%32s = %", n_digits,"d\n"),
                "Number of standalone items",
                n_standalone_items))
      text_after_df <- paste0(
        text_after_df,
        sprintf(paste0("%32s = %", n_digits,"d\n"),
                "Total number of items", n_all_items))
    } else {
      text_after_df <- paste0("Total number of items = ", length(x), "\n")
    }
    return(text_after_df)
  }
  if (has_testlet) result$testlet[is.na(result$testlet)] <- ""

  ## Print first n rows ##
  # if number of items are between 1-20, print all items. If number of items
  # are larger than 20, print first 10 items.
  n_all_items <- x$n$items
  if (is.null(n)) n <- ifelse(n_all_items <= 20, n_all_items, 10)

  if (n < 1) {
    result <- ""
    text_after_df <- print_ip_size(x)
  } else if (n < x$n$items) {
    result <- result[1:n, ]
    text_after_df <- paste0(
      format_text(paste0("# ... with ", x$n$items - n, " more items\n"),
                  fg = "light gray", italic = TRUE),
      print_ip_size(x))
  } else {
    text_after_df <- paste0()
  }

  print(result)
  cat(text_after_df)
}

###############################################################################@
############################# print.Itempool ###################################
###############################################################################@
#' Show an \code{\link{Itempool-class}} object
#'
#' @param x An \code{\link{Itempool-class}} object that will be showed.
#' @param ... Additional arguments. For example, an argument \code{n = 14},
#'   will print 14 items to the console.
#'
#' @export
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
#'
#' @examples
#' ip <- generate_ip(n = 42)
#' print(ip)
#' print(ip, n = 3)
#' print(ip, n = 12)
#' print(ip, n = Inf)
#'
setMethod("print", "Itempool", function(x, ...)  {
  args <- list(...)
  .print.Itempool(x = x, n = switch("n" %in% names(args), args$n, NULL))
  })


###############################################################################@
############################# show.Itempool ###################################
###############################################################################@
#' Show an \code{\link{Itempool-class}} object
#'
#' @param object An \code{\link{Itempool-class}} object that will be showed.
#'
#' @export
#'
#' @rdname show
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
setMethod("show", "Itempool", function(object) {.print.Itempool(object)})


###############################################################################@
############################# as.Itempool #####################################
###############################################################################@
#' Coerce a given object to \code{\link{Itempool-class}} object
#'
#'
#' @description This function is a wrapper for \code{\link{itempool}} function.
#' It is recommended to use that function.
#'
#' @param ... The object that is desired to be converted to an  'Itempool'
#'          object. Also additional arguments related to the Itempool.
#'
#' @return An \code{\link{Itempool-class}} object.
#'
#' @include Item-class.R
#' @include Item-class-methods.R
#' @include Itempool-class.R
#'
#'
#' @seealso \code{\link{itempool}}
#'
#' @export
#' @importFrom methods validObject
#'
#' @author Emre Gonulates
"as.Itempool" <- function(...)
{
  args <- list(...)
  return(do.call("itempool", args))
}


###############################################################################@
############################# convert_model (generic)  #########################
###############################################################################@
#' Convert model parameters from one model to another
#'
#' @name convert_model
#' @description
#' This is especially handy for converting IRT models with less parameters
#' (such as 1 parameter logistic model) to higher dimensional models such as
#' three parameter logistic model.
#'
#' @param ip An \code{\link{Item-class}} or \code{\link{Itempool-class}} object
#' @param target_model The target model that the conversion will be made.
#'
#' @include Item-class.R
#'
#' @rdname convert_model
#'
#' @author Emre Gonulates
#'
#' @return An 'Item' object with new model parameters will be returned.
setGeneric("convert_model",
           def = function(ip, target_model = "3PL")
           {standardGeneric ("convert_model")},
           useAsDefault =
             function(ip, target_model = "3PL")
             stop(paste0("Invalid object type. Only 'Item', 'Testlet' or ",
                          "'Itempool' types objects can be converted. "))
)


###############################################################################@
############################# convert_model (Item) #############################
###############################################################################@
#' Convert model parameters from one model to another.
#'
#' @rdname convert_model
#'
#' @export
#'
#' @importFrom methods new
#'
#' @author Emre Gonulates
#'
setMethod(
  f = "convert_model", signature = c(ip = "Item"),
  function(ip, target_model = "3PL") {
    uirt_models <- names(Pmodels)[sapply(Pmodels, function(x)
      x$model_family == "UIRT")]
    available_target_models <- c(uirt_models, "GPCM")
    if (!target_model %in% available_target_models)
      stop(paste0("Invalid target_model. Target model should be either: ",
                  paste0("'", available_target_models, "'", collapse = ", "),
                  "."))
    current_model <- ip$model
    current_pars <- ip$parameters
    # Both target and item belongs to unidimensional irt models
    if (target_model %in% uirt_models && current_model %in% uirt_models) {
      switch(
        target_model,
        "Rasch" = {
          output <- new("Item", model = "Rasch",
                        parameters = list(b = 0),
                        se_parameters = ip@se_parameters, id = ip@id,
                        content = ip@content,
                        misc = ip@misc
                        )
        },
        "1PL" = {
          output <- new("Item", model = "1PL",
                        parameters = list(b = 0, D = default_D_scaling),
                        se_parameters = ip@se_parameters, id = ip@id,
                        content = ip@content,
                        misc = ip@misc
                        )
        },
        "2PL" = {
          output <- new(
            "Item", model = "2PL",
            parameters = list(a = 1, b = 0, D = default_D_scaling),
            se_parameters = ip@se_parameters, id = ip@id, content = ip@content,
            misc = ip@misc)
        },
        "3PL" = {
          output <- new(
            "Item", model = "3PL",
            parameters = list(a = 1, b = 0, c = 0, D = default_D_scaling),
            se_parameters = ip@se_parameters, id = ip@id, content = ip@content,
            misc = ip@misc)
        },
        "4PL" = {
          output <- new(
            "Item", model = "4PL",
            parameters = list(a = 1, b = 0, c = 0, d = 1, D = default_D_scaling),
            se_parameters = ip@se_parameters, id = ip@id, content = ip@content,
            misc = ip@misc)
        })
      # Create a dummy Item
      overlappingPars <- names(current_pars)[
        names(current_pars) %in% names(output@parameters)]
      output@parameters[overlappingPars] <- current_pars[overlappingPars]
    } else if (current_model == "GPCM2" && target_model == "GPCM") {
      output <- new(
        "Item", model = "GPCM",
        parameters = list(a = current_pars$a, b = current_pars$b-current_pars$d,
                          D = current_pars$D),
        se_parameters = ip@se_parameters, id = ip@id, content = ip@content,
        misc = ip@misc)

    } else
      stop("The Item cannot be converted to the target model.")
    return(output)
  })

###############################################################################@
############################# convert_model (Itempool) ########################
###############################################################################@
#' Convert model parameters from one model to another.
#'
#' @rdname convert_model
#'
#' @keywords internal
#'
#' @export
#'
#' @importFrom methods new
#'
#' @author Emre Gonulates
#'
setMethod(
  f = "convert_model", signature = c(ip = "Itempool"),
  function(ip, target_model = "3PL")
  {
    return(new(
      Class = "Itempool", item_list = lapply(
        ip@item_list, FUN = function(y) convert_model(
          ip = y, target_model = target_model))))
  })


###############################################################################@
############################# convert_model (Testlet) ##########################
###############################################################################@
#' Convert model parameters from one model to another.
#'
#' @rdname convert_model
#'
#' @keywords internal
#'
#' @export
#'
#' @author Emre Gonulates
#'
setMethod(
  f = "convert_model", signature = c(ip = "Testlet"),
  function(ip, target_model = "3PL")
  {
    ip@item_list <- convert_model(ip = ip@item_list,
                                  target_model = target_model)
    return(ip)
  })


###############################################################################@
############################# generate_item ####################################
###############################################################################@
#' Generate a random Item object
#'
#' @param model The model of the Item object.
#' @param n_categories For polytomous models, the number of categories for an
#'   'item' object.
#' @param se_parameters The values of parameter standard errors, i.e. a list
#'   object with elements named as parameter names (excluding \code{"D"}
#'   parameter).
#'
#'   If the value is \code{TRUE}, this function will generate standard error
#'   values from a uniform distribution between 0.05 and 0.75 for each parameter.
#' @param ... Additional parameters passed to \code{itempool()} function.
#'
#' @return An \code{\link{Item-class}} object
#'
#' @importFrom stats rnorm rlnorm runif
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @examples
#' # By default, a '3PL' model Item generated
#' generate_item()
#' # Generate item pools for other models
#' generate_item("Rasch")
#' generate_item("1PL")
#' generate_item("2PL")
#' generate_item("4PL")
#' # Polytomous items
#' generate_item("GRM")
#' generate_item("GPCM")
#' generate_item("PCM")
#' generate_item("GPCM2")
#' # Different number of categories
#' generate_item("GRM", n_categories = 2)
#' generate_item("GPCM", n_categories = 5)
#'
#' # Generate standard errors for item parameters
#' generate_item(se_parameters = TRUE)
#'
generate_item <- function(model = "3PL", n_categories = 4,
                          se_parameters = NULL, ...) {
  ##### Argument Checks #@###
  if (!is.numeric(n_categories) || length(n_categories) != 1 ||
      n_categories < 2) {
    stop("Invalid number of categories, 'n_categories'. Provide an integer ",
         "larger than 1.")
  } else n_categories <- as.integer(n_categories)

  if (is.null(model) || length(model) != 1 || !is.character(model) ||
      !model %in% names(Pmodels))
    stop("Invalid model argument.")

  pars <- list()
  if (model %in%
      names(Pmodels)[sapply(Pmodels, function(x) x$model_family == "UIRT")]) {
    pars$b <- round(rnorm(1), 4)
    if (model %in% c("1PL", "2PL", "3PL", "4PL")) {
      pars$D <- default_D_scaling
      if (model %in% c("2PL", "3PL", "4PL")) {
        pars$a <- round(rlnorm(1, 0, 0.3), 4)
        if (model %in% c("3PL", "4PL")) {
          pars$c <- round(runif(1, 0, 0.3), 4)
          if (model == "4PL") pars$d <- round(runif(1, 0.90, 0.99), 4)
        }
      }
    }
  } else if (model %in%
      names(Pmodels)[sapply(Pmodels, function(x) x$model_family == "PIRT")]) {
    # Create threshold parameters
    while (TRUE) {
      pars$b <- sort(round(rnorm(n_categories - 1, sd = .75), 4))
      temp_diff <- pars$b[-1] - pars$b[-(n_categories - 1)]
      if (n_categories == 2 || min(temp_diff) > 0.3)
        break
    }
    # For GPCM2, switch b to d and create a new b location parameter
    if (model == "GPCM2") {
      pars$d <- pars$b
      pars$b <- rnorm(1, 0, 0.33)
    }

    if (model %in% c("GRM", "GPCM", "GPCM2")) {
      pars$a <- round(rlnorm(1, 0, 0.3), 4)
      pars$D <- default_D_scaling
    }
  } else {
    warning("This model has not been implemented yet.")
    return(NULL)
  }

  # Create se_parameters if se_parameters value is `TRUE`
  if (!is.null(se_parameters) && is.logical(se_parameters) &&
      length(se_parameters) == 1 && se_parameters) {
    # se_parameters <- stats::setNames(
    #   vector("list", length(Pmodels[[model]]$se_parameters)),
    #   Pmodels[[model]]$se_parameters)
    se_parameters <- names(Pmodels[[model]]$parameters)[
          sapply(Pmodels[[model]]$parameters, `[[`, "se")];
    se_parameters <- stats::setNames(
      vector("list", length(se_parameters)), se_parameters)
    for (par_name in names(se_parameters))
      se_parameters[[par_name]] <- runif(length(pars[[par_name]]), min = 0.05,
                                         max = 0.75)
  }

  return(item(parameters = pars, model = model, se_parameters = se_parameters,
              ...))
}


###############################################################################@
############################# generate_ip ######################################
###############################################################################@
#' Generate a random \code{Itempool} object
#'
#' @param model The model of the item pool
#' @param n The number of items in the item pool.
#' @param output The type of object returned. The default value is
#'   \code{"Itempool"}.
#'   \describe{
#'     \item{\strong{\code{"Itempool"}}}{Return an
#'       \code{\link{Itempool-class}} object.
#'       }
#'     \item{\strong{\code{"Item"}}}{If \code{n = 1} return an
#'       \code{\link{Item-class}} object. If \code{n > 1}, returns a list of
#'       \code{\link{Item-class}} object.
#'       }
#'     \item{\strong{\code{"list"}}}{Return a list of item
#'       \code{\link{Item-class}} objects.
#'       }
#'     }
#' @param n_categories For polytomous items, designate the number of categories
#'   each item should have. It can be a single integer value larger than 1. In
#'   this case all of the polytomous items will have this number of categories.
#'   It can be a vector of length \code{n} designating the categories of each
#'   item. For dichotomous items, the values in \code{n_categories} will be
#'   ignored.
#' @param se_parameters The values of parameter standard errors for each item,
#'   i.e. a list
#'   object with elements named as parameter names (excluding \code{"D"}
#'   parameter).
#'
#'   If the value is \code{TRUE}, this function will generate standard error
#'   values from a uniform distribution between 0.05 and 0.75 for each
#'   parameter of each item.
#' @param ... Additional parameters passed to \code{itempool()} function.
#'
#' @return An \code{\link{Itempool-class}} object
#'
#' @author Emre Gonulates
#'
#' @export
#'
#' @examples
#' # By default, a '3PL' model item pool generated
#' generate_ip()
#' # Designate the number of items
#' generate_ip(n = 12)
#' # Generate item pools for other models
#' generate_ip(model = "Rasch")
#' generate_ip(model = "1PL")
#' generate_ip(model = "2PL")
#' generate_ip(model = "4PL")
#' generate_ip(model = "GRM") # Graded Response Model
#' generate_ip(model = "GPCM") # Generalized Partial Credit Model
#' generate_ip(model = "PCM") # Partial Credit Model
#' generate_ip(model = "GPCM2") # Reparametrized GPCM
#' # Mixture of models
#' generate_ip(model = c("4PL", "Rasch"))
#' generate_ip(model = sample(c("4PL", "GPCM"), 12, TRUE))
#' generate_ip(model = c("2PL", "GRM", "Rasch"), n = 11)
#'
#' # Generate parameters standard errors for each item
#' generate_ip(se_paramters = TRUE)
#'
#' # Generate an item pool consist of testlets and standalone items
#' temp_list <- list(ids = paste0("testlet-", 1:7), n = c(2, 3, 4, 2, 3, 4, 2))
#' ip <- itempool(sample(c(generate_ip(n = 10, output = "list"),
#'                         sapply(1:length(temp_list$id), function(i)
#'                           generate_testlet(id = temp_list$id[i],
#'                                            n = temp_list$item_models[i])))))
#'
generate_ip <- function(model = "3PL", n = NULL, output = "Itempool",
                        n_categories = 4, se_parameters = NULL, ...) {
  if (is.null(n)) n <- ifelse(length(model) > 1, length(model),
                              sample(10:20, 1))
  model <- rep(model, length.out = n)

  if (is.null(n_categories) || !length(n_categories) %in% c(1, n))
    stop("Invalid n_categories value.")
  n_categories <- rep(n_categories, length.out = n)

  item_list <- lapply(1:n, function(i) generate_item(
    model = model[i], n_categories = n_categories[i],
    # The following returns TRUE if se_parameters and NULL otherwise
    se_parameters = switch(is.list(se_parameters)+1, se_parameters, NULL))
    )
  # If se_parameters is TRUE, no need to add se_parameters to itempool()
  # function because the se_parameters are already set in item_list.
  if (!is.null(se_parameters) && is.logical(se_parameters) &&
      length(se_parameters) == 1 && se_parameters) {
    ip <- itempool(item_list, ...)
  } else {
    ip <- itempool(item_list, se_parameters = se_parameters, ...)
  }

  if (output == "Itempool") {
    return(ip)
  } else if (output == "list" || (output == "Item" && n > 1)) {
    return(ip@item_list)
  } else if (output == "Item") {
    return(ip@item_list[[1]])
  } else stop("Invalid 'output' value.")
}

###############################################################################@
############################# generate_testlet #################################
###############################################################################@
#' Generate a random Testlet object
#'
#' @param model The model of the Testlet
#' @param n The number of items in the Testlet.
#' @param item_models A single model name or a vector of model names with the
#'   size of n that represents the models of items in the Testlet object.
#' @param item_id_preamble The preamble for the item ids within the Testlet.
#' @param n_categories For polytomous items, designate the number of categories
#'   each item should have. It can be a single integer value larger than 1. In
#'   this case all of the polytomous items of the testlet will have this number
#'   of categories. It can be a vector of length \code{n} designating the
#'   categories of each item. For dichotomous items, the values in
#'   \code{n_categories} will be ignored.
#' @param ... Additional parameters passed to \code{testlet()} function.
#'
#' @return A \code{\link{Testlet-class}} object
#'
#' @author Emre Gonulates
#'
#' @export
#'
#' @examples
#' # By default, a Testlet object with '3PL' model items generated
#' generate_testlet()
#' # Designate the number of items in the testlet
#' generate_testlet(n = 12)
#' # Set the id of the testlet
#' generate_testlet(id = "my-testlet")
#' # Designate the id of testlet and preamble for item ids
#' generate_testlet(id = "my-testlet", item_id_preamble = "mt-")
#' # Generate item pools for other models
#' generate_testlet(item_model = "Rasch")
#' generate_testlet(item_model = "1PL")
#' generate_testlet(item_model = "2PL")
#' generate_testlet(item_model = "4PL")
#' generate_testlet(item_model = "GRM") # Graded Response Model
#' generate_testlet(item_model = "GPCM") # Generalized Partial Credit Model
#' generate_testlet(item_model = "PCM") # Partial Credit Model
#' generate_testlet(item_model = "GPCM2") # Reparametrized GPCM
#' # Mixture of models
#' generate_testlet(item_models = c("4PL", "Rasch"))
#' generate_testlet(model = c("2PL", "GRM", "Rasch"), n = 11)
#'
#' # Generating multiple testlet objects with custom ids
#' sapply(paste0("testlet-", 1:4), function(x) generate_testlet(id = x))
#'
#' # Generate testlet with dichotomous and polytomous with different number of
#' # categories.
#' generate_testlet(
#'   item_models = c("3PL", "GRM", "GPCM", "GRM", "2PL"),
#'   n_categories = c(2, 3, 6, 7, 2))
#'
#' # # Generating multiple testlet objects with custom ids and item models and
#' # # put them in an item pool:
#' # temp_list <- list(ids = paste0("testlet-", 1:3),
#' #                   item_models = c("Rasch", "2PL", "GPCM"))
#' # itempool(sapply(1:length(temp_list$id), function(i)
#' #   generate_testlet(id = temp_list$id[i],
#' #   item_models = temp_list$item_models[i])))
#'
generate_testlet <- function(model = "BTM", n = NULL,
                             item_models = "3PL", item_id_preamble = NULL,
                             n_categories = 4, ...) {
  args <- list(...)
  if (is.null(n)) {
    n <- sample(2:5, 1)
    if (length(item_models) > 1) n <- length(item_models)
  }
  if (is.null(item_models) || !is.character(item_models)) {
    item_models <- rep("3PL", length.out = n)
  } else if (length(item_models) == 1 || length(item_models) != n) {
    item_models <- rep(item_models, length.out = n)
  }
  ip <- generate_ip(model = item_models, n_categories = n_categories)
  # item_list <- lapply(item_models, generate_item)
  # ip <- itempool(item_list)
  if (is.null(item_id_preamble)) {
    if ("id" %in% names(args)) {
      ip$id <- paste0(args$id, "-", ip$id)
    } else
      ip$id <- paste0("t-", ip$id)
  } else {
      ip$id <- paste0(item_id_preamble, ip$id)
  }
  return(testlet(ip, ...))
}



###############################################################################@
############################# add_misc (generic)  ##############################
###############################################################################@
#' Add or change a named value to 'misc' slot of an \code{\link{Item-class}},
#' \code{\link{Itempool-class}} or \code{\link{Testlet-class}} object.
#'
#' @name add_misc
#'
#' @param ip An \code{\link{Item-class}}, \code{\link{Testlet-class}} or
#'   \code{\link{Itempool-class}} object.
#' @param value A list where each element should be named. Elements within the
#'   list will be added to 'misc' slot.
#'
#' @return An object with added 'misc' slot.
#'
#' @include Item-class.R
#'
#' @rdname add_misc
#'
#' @author Emre Gonulates
#'
#' @export
#'
#' @examples
#' item <- item(b = 1)
#' add_misc(item, list(sympson_hetter_k = .75))
#'
setGeneric("add_misc",
           def = function(ip, value)
           {standardGeneric ("add_misc")},
           useAsDefault =
             function(ip, value)
             stop(paste0("Invalid object type. Only 'Item', 'Testlet' type or ",
                          "'Itempool' type objects can be used. "))
)


# The following is a generic function of 'add_misc' since it is basically the
# same for 'Item', 'Itempool' or 'Testlet'
.add_misc <- function(ip, value) {
  if (!is(value, "list") || is.null(names(value)) || any(names(value) == ""))
    stop("The 'value' argument should be a named list element.")
  if (is.null(ip@misc)) ip@misc = list()
  for (i in seq_along(length(value)))
    ip@misc[[names(value[i])]] <- value[[i]]
  return(ip)
}

###############################################################################@
############################# add_misc (Item) ##################################
###############################################################################@
#' @rdname add_misc
#'
#' @export
#'
setMethod(
  f = "add_misc", signature = c(ip = "Item"),
  function(ip, value) {
    return(.add_misc(ip, value))
  })

###############################################################################@
############################# add_misc (Testlet) ###############################
###############################################################################@
#' @rdname add_misc
#'
#' @export
#'
setMethod(
  f = "add_misc", signature = c(ip = "Testlet"),
  function(ip, value) {
    return(.add_misc(ip, value))
  })


###############################################################################@
############################# add_misc (Itempool) #############################
###############################################################################@
#' @rdname add_misc
#'
#' @export
#'
setMethod(
  f = "add_misc", signature = c(ip = "Itempool"),
  function(ip, value) {
    return(.add_misc(ip, value))
  })

