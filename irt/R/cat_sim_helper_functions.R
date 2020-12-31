
###############################################################################@
############################# print.cat_design ##################################
###############################################################################@
#' Prints cat_design objects.
#'
#' @param x A \code{cat_design} object.
#' @param ... further arguments passed to or from other methods.
#' @param verbose If \code{TRUE}, a list object will be returned for each
#'          step.
#'
#' @keywords internal
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @method print cat_design
#'
#' @examples
#' ip <- generate_ip(n = 5)
#' cd <- create_cat_design(ip = ip, next_item_rule = 'random',
#'                         termination_rule = 'min_item',
#'                         termination_par = list('min_item' = 5))
#' cd
print.cat_design <- function(x, ..., verbose = FALSE) {
  if (verbose) {
    NextMethod()
    return(NULL)
  }
  # This function extracts the sub parameters from 'ability_est_par' or
  # 'termination_par'.
  # 'g_par_name': is the general parameter name. It should be either
  #    'ability_est_par'  or 'termination_par'
  # 'title': The title that will be printed for common section.
  # 'preamble_for_df': this is the preamble that will be put before the
  #   parameter name in steps_df data frame.
  process_parameters <- function(g_par_name, steps_df, steps_list,
                                 title, preamble_for_df) {
    # Extract the prior parameter names for each step
    par_names <- lapply(x$step, function(y) names(y[[g_par_name]]))
    par_names_union <- unique(unlist(par_names))
    par_names_intersection <- Reduce(intersect, par_names)
    common_parameters_exists <- FALSE
    if (length(par_names_intersection) > 0)
      for (pp_name in par_names_intersection)
        # Extract all parameters with the name pp_name
        # Check whether all of parameters are the same
        if (length(unique(lapply(x$step,
                                 function(y) y[[g_par_name]][[pp_name]]))) == 1)
          common_parameters_exists <- TRUE
    if (common_parameters_exists)
      cat(paste0(title, ": \n"))
    if (!is.null(par_names_intersection))
      max_char <- max(nchar(par_names_intersection))
    # Deal with the common prior parameter elements, either put them in the
    # common parameter section or put them in steps_df
    for (pp_name in par_names_intersection) {
      # Extract all parameters with the name pp_name
      temp <- lapply(x$step, function(y) y[[g_par_name]][[pp_name]])
      # Check whether all of parameters are the same
      if (length(unique(temp)) == 1) {
        temp <- unique(temp)[[1]]
        if (class(temp) %in% printable_classes) {
          cat(paste0("    ", sprintf(paste0("%", max_char, "s"), pp_name), ": ",
                     print_with_quotes(temp), "\n"))
          # cat(paste0("    ", pp_name, ": ", print_with_quotes(temp), "\n"))
        } else if (is(temp, "list") &&
                   all(sapply(temp, class) %in% printable_classes)) {
          cat(sprintf(paste0("    %", max_char, "s: \n"), pp_name))
          # cat(sprintf("    %s: \n", pp_name))
          for (i in seq_len(length(temp)))
            cat(paste0("        ", sprintf(paste0("%", max_char, "s:"),
                                           names(temp)[i]), ": ",
                       print_with_quotes(temp[[i]]), "\n"))
            # cat(paste0("        ", names(temp)[i], ": ",
            #            print_with_quotes(temp[[i]]), "\n"))
        } else {
          cat(sprintf(paste0("    %", max_char, "s: '\n"), pp_name))
          # cat(sprintf("    %s: '\n", pp_name))
          # cat(sprintf("    %s: '\n", pp_name))
          print(temp)
        }
      } else {
        # If all of them are printable classes add them to the steps_df otherwise
        # add them to steps_list
        if (all(sapply(temp, class) %in% printable_classes)) {
          steps_df[, paste0(preamble_for_df, "_", pp_name)] <- unlist(temp)
        } else {
          for (i in seq_len(number_of_steps))
            steps_list[[i]][pp_name] <- temp[[i]]
        }
      }
    }
    # Deal with the elements that are not common accross the steps
    for (pp_name in setdiff(par_names_union, par_names_intersection)) {
      # Add them directly to steps_list
      for (i in seq_len(number_of_steps))
        steps_list[[i]][pp_name] <- x$step[[i]][[g_par_name]][pp_name]
    }
    return(list(steps_df = steps_df, steps_list = steps_list))
  }

  print_with_quotes <- function(x)
    ifelse(test = is(x, "character"), yes = paste0("'", x, "'"), no = paste0(x))

  number_of_steps <- length(x$step)
  steps_df <- data.frame(step = seq_len(number_of_steps))
  steps_list <- vector("list", number_of_steps)
  printable_classes <-  c('logical', 'character', 'integer', 'numeric')

  cat("CAT Design\n")
  cat("--------------------------------------------------\n")
  if (!is.null(x$title))
    cat(sprintf("Title: '%s'\n", x$title))
  cat(sprintf("Item Pool Size: %d\n", length(x$ip)))
  cat(sprintf("Maximum Test Length: %d\n", x$max_test_length))

  ## First Item Parameters ##
  cat(sprintf("First Item Rule: '%s'\n", x$first_item_rule))
  cat("First Item Parameters: \n")
  for (i in seq_len(length(x$first_item_par)))
    if (class(x$first_item_par[[i]]) %in% printable_classes) {
      cat(paste0("   ", names(x$first_item_par)[i], ": ",
                 print_with_quotes(x$first_item_par[[i]]), "\n"))
    } else {
      cat(sprintf("    %s: \n", names(x$first_item_par)[i]))
      print(x$first_item_par[[i]])
    }

  # Next item rule
  temp <- vapply(x$step, function(y) y$next_item_rule, character(1))
  if (all(temp == temp[1])) {
    cat(sprintf("Next Item Rule: %s\n", print_with_quotes(temp[1])))
  } else steps_df$next_item_rule <- temp

  temp <- process_parameters(g_par_name = 'next_item_par',
                             steps_df = steps_df,
                             steps_list = steps_list,
                             title = "Next Item Parameters",
                             preamble_for_df = "NIP")
  steps_df <- temp$steps_df
  steps_list <- temp$steps_list

  # Ability Estimation Rule
  temp <- vapply(x$step, function(y) y$ability_est_rule, character(1))
  if (all(temp == temp[1])) {
    cat(sprintf("Ability Estimation Rule: %s\n", print_with_quotes(temp[1])))
  } else steps_df$ability_est_rule <- temp

  ## Ability Estimation Parameters ##
  temp <- process_parameters(g_par_name = 'ability_est_par',
                             steps_df = steps_df, steps_list = steps_list,
                             title = "Ability Estimation Parameters",
                             preamble_for_df = "AEP")
  steps_df <- temp$steps_df
  steps_list <- temp$steps_list

  ## Final Ability Estimation Parameters ##
  if (!is.null(x$final_ability_est_rule)) {
    cat(sprintf("Final Ability Estimation Rule: '%s'\n", x$final_ability_est_rule))
    cat("Final Ability Estimation Parameters: \n")
    for (i in seq_len(length(x$final_ability_est_par)))
      if (class(x$final_ability_est_par[[i]]) %in% printable_classes) {
        cat(paste0("    ", names(x$final_ability_est_par)[i], ": ",
                   print_with_quotes(x$final_ability_est_par[[i]]), "\n"))
      } else {
        cat(sprintf("    %s: \n", names(x$final_ability_est_par)[i]))
        print(x$final_ability_est_par[[i]])
      }
  }

  ## Test Termination Rule ##
  cat("Test Termination Rules and Parameters (in order): \n")
  max_char <- max(nchar(unlist(sapply(x$termination_par, names))))
  for (i in seq_len(length(x$termination_par)))
    if (length(x$termination_par[[i]]) > 1) {
      cat(paste0("    (Rule ", i, ") ", names(x$termination_par)[i], ": \n"))
      for (j in 1:length(x$termination_par[[i]])) {
        cat(paste0("        ", sprintf(paste0("%", max_char, "s: "),
                                       names(x$termination_par[[i]])[j]),
                   print_with_quotes(x$termination_par[[i]][[j]]), "\n"))
      }
    } else if (length(x$termination_par[[i]]) == 1) {
      cat(paste0("    (Rule ", i, ") ", names(x$termination_par)[i], ": ",
                 print_with_quotes(x$termination_par[[i]]),
                 "\n"))

    }

    # if (class(x$termination_par[[i]]) %in% printable_classes) {
    #   cat(paste0("   (Rule ", i, ") ", names(x$termination_par)[i], ": ",
    #              print_with_quotes(x$termination_par[[i]]), "\n"))
    # } else {
    #   cat(sprintf("    (Rule %d) %s: \n", i, names(x$termination_par)[i]))
    #   print(x$termination_par[[i]])
    # }

  # Print irregular parameters
  if (ncol(steps_df) > 1) {
    cat("\nStep Arguments that are not Common:\n")
    print(steps_df)
  }
  for (i in number_of_steps:1)
    if (is.null(steps_list[[i]])) steps_list[[i]] <- NULL
  if (length(steps_list) > 0) {
    cat("\nUncommon and Irregular Step Arguments:\n")
    print(steps_list)
  }
}

###############################################################################@
############################# summary.cat_output ###############################
###############################################################################@
#' Summarizes the raw output of cat_sim
#'
#' @description This function summarizes a list consist of cat_output objects.
#' It returns a summary data frame of the CAT simulation.
#'
#' @param object This is a cat_output object or a list object containing
#'   elements that are "cat_output" class.
#' @param ... Additional arguments.
#' @param cols The variables that will be included in the summary. There should
#'          be at least one column. Available columns are:
#'          \describe{
#'            \item{true_ability}{True ability of the simulee}
#'            \item{est_ability}{Ability Estimate}
#'            \item{se}{Standard Error of the ability estimate}
#'            \item{test_length}{Test length.}
#'            \item{bias}{The difference between true ability and ability
#'              estimate}
#'            \item{mse}{Mean squared error}
#'          }
#' @return This function returns a summary data frame of adaptive tests. Each
#' row will represent a different adaptive test.
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @seealso \code{\link{cat_sim}}
#'
#' @examples
#' n <- 100 # number of items
#' ip <- generate_ip(n = n,
#'                   content = sample(c("Algebra", "Arithmetic", "Geometry"),
#'                                    n, replace = TRUE))
#' cd <- create_cat_design(ip = ip, next_item_rule = 'mfi',
#'                         termination_rule = 'max_item',
#'                         termination_par = list(max_item = 10))
#' cat_data <- cat_sim(true_ability = rnorm(5), cd = cd)
#' summary(cat_data)
summary.cat_output <- function(
  object, ...,
  cols = c("true_ability", "est_ability", "se", "test_length")) {
  # Check whether it is one or multiple CAT outputs
  if (is(object, "cat_output")) { # There is a single cat_output element
     nrows <- 1
  } else if (all(sapply(object, is, "cat_output"))) { # a list of cat_output
    nrows <- length(object)
  } else
    stop("All of the elements of 'cat_output' should be 'cat_output' class.")
  ncols <- length(cols)

  # if (!all(sapply(x, FUN = function(y) is(y, "cat_output"))))
  #   stop("All of the elements of 'cat_output' should be 'cat_output' class.")
  if (length(cols) < 1) stop("There should be at least one column.")
  if (!all(cols %in% c("true_ability", "est_ability", "se", "test_length",
                       "bias", "mse")))
    stop("Inadmissable column. 'cols' should be composed of one of the
         following: 'true_ability', 'est_ability', 'se', 'test_length', 'bias',
         'mse'.")
  cat_summary <- data.frame(matrix(vector(), nrows, ncols,
                                   dimnames = list(c(), cols)),
                            stringsAsFactors = FALSE)
  # Function returns the true ability
  get_true_ability <- function() {
    if (nrows == 1) {
      return(object$true_ability)
      } else
      return(unlist(sapply(object, `[[`, "true_ability")))
  }
  get_est_ability <- function() {
    if (nrows == 1) {
      eh <- object$est_history
      return(eh[[length(eh)]]$est_after)
    } else {
      return(sapply(object, function(x) {
        eh <- x$est_history; eh[[length(eh)]]$est_after}))
    }
  }
  for (col in cols) {
    switch(
      col,
      "true_ability" = {
        cat_summary$true_ability <- get_true_ability()
      },
      "est_ability" = {
        cat_summary$est_ability <- get_est_ability()
      },
      "se" = {
        if (nrows == 1) {
          eh <- object$est_history
          cat_summary$se <- eh[[length(eh)]]$se_after
        } else {
          cat_summary$se <- sapply(
            object, function(x) {
              eh = x$est_history; eh[[length(eh)]]$se_after
              })
        }
      },
      "test_length" = {
        if (nrows == 1) {
          cat_summary$test_length <- length(object$est_history)
        } else {
          cat_summary$test_length <- sapply(
            object, function(x) {length(x$est_history)})
        }
      },
      "bias" = {
        cat_summary$bias <- get_est_ability() - get_true_ability()
      },
      "mse" = {
        cat_summary$mse <- (get_est_ability() - get_true_ability())^2
        }
    )
  }
  return(cat_summary)
}


###############################################################################@
############################# .print.cat_output ################################
###############################################################################@
#' Prints the raw output of cat_sim
#' @description This function prints a data frame that shows all of the steps of
#'   a CAT for a single examinee.
#'
#' @param x This is a cat_output object which has "cat_output" class.
#' @param ... Additional arguments.
#' @param silent If TRUE, no output will be printed on the console, only
#'   a data frame will be returned.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @author Emre Gonulates
#'
#' @seealso \code{\link{cat_sim}}
#'
.print.cat_output <- function(x, ..., silent = FALSE) {
  # Check whether it is one or multiple CAT outputs
  if (!is(x, "cat_output"))
    stop("x should be 'cat_output' class.")
  est_history <- x$est_history
  n_items <- length(est_history)

  output <- data.frame(
    est_before = sapply(est_history, `[[`, "est_before"),
    se_before = sapply(est_history, `[[`, "se_before"),
    testlet_id = sapply(est_history, function(i)
      ifelse(is.null(i$testlet), NA, i$testlet$id)),
    item_id = sapply(est_history, function(i)
      ifelse(is.null(i$item), NA, i$item$id)),
    resp = sapply(est_history, `[[`, "resp"),
    est_after = sapply(est_history, `[[`, "est_after"),
    se_after = sapply(est_history, `[[`, "se_after"),
    stringsAsFactors = FALSE
  )

  if (!silent) {
    cat("An object of class 'cat_output'.\n")
    cat(paste0("True Ability: ", x$true_ability, "\n\n"))
    print(output)
  }
  invisible(output)
}

###############################################################################@
############################# show.cat_output ##################################
###############################################################################@
#' This method shows an "cat_output" class object
#'
#' @param object An 'cat_output' class object that will be showed.
#'
#' @export
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
show.cat_output <- function(object) .print.cat_output(object)


###############################################################################@
############################# print.cat_output #################################
###############################################################################@
#' This method prints an "cat_output" class object
#'
#' @param object An 'cat_output' class object that will be printed.
#' @param ... Additional arguments
#' @param silent If TRUE, no output will be printed on the console, only
#'   a data frame will be returned.
#'
#' @export
#'
#' @method print cat_output
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
print.cat_output <- function(x, ..., silent = FALSE)
  .print.cat_output(x, ..., silent = silent)


###############################################################################@
############################# summary.list #####################################
###############################################################################@
#' If a list object consists of all "cat_output" objects, then it will run
#' summary.cat_output.
#'
#' @param object A list object consists of all "cat_output" objects.
#' @param ... Arguments passed to the \code{summary.cat_output()} function.
#' @return A data frame that summarizes the CAT outputs.
#'
#' @export
#' @keywords internal
#'
#' @author Emre Gonulates
#'
summary.list <- function(object, ...) {
  if (all(sapply(object, class) == "cat_output")) {
    summary.cat_output(object, ...) } else NextMethod()
}


###############################################################################@
############################# get_cat_response_data ############################
###############################################################################@
#' Extracts the response data of CAT output.
#'
#' @description This function extracts the response data from a single
#' \code{cat_output} object or a list of \code{cat_output} objects and gives
#' either a vector (if there is a single \code{cat_output} object) or a matrix
#' (if there is a lit of \code{cat_output} objects) of response data.
#'
#' If \code{cd}, cat design, object is given, then the item pool in the
#' \code{cd} will be used.
#'
#'
#' @param cat_sim_output This is a list object containing elements that are
#' "cat_output" class.
#' @param cd A \code{cat_design} object that is created by function
#'          \code{create_cat_design}.
#' @param remove_na If \code{TRUE}, the columns that are all \code{NA} will be
#'          removed.
#' @param attach_summary If \code{TRUE}, the summary of each CAT will be
#'   attached to the beginning of the response string as columns. The default
#'   value is \code{FALSE}.
#' @return This function returns a response matrix of adaptive tests. If the
#' input is a list of \code{cat_output}, then the rows will represent examinees
#' and columns will represent items. For single \code{cat_output} object the
#' vector names will be the element
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @seealso \code{\link{cat_sim}}
#'
#' @examples
#' n <- 40 # number of items
#' ip <- generate_ip(n = n,
#'                   content = sample(c("Algebra", "Arithmetic", "Geometry"),
#'                                    n, replace = TRUE))
#' cd <- create_cat_design(ip = ip, next_item_rule = 'mfi',
#'                         termination_rule = 'max_item',
#'                         termination_par = list(max_item = 10))
#' cat_data <- cat_sim(true_ability = rnorm(10), cd = cd)
#' get_cat_response_data(cat_sim_output = cat_data, cd)
#'
get_cat_response_data <- function(cat_sim_output, cd = NULL,
                                  remove_na = FALSE, attach_summary = FALSE) {
  # Check whether it is one or multiple CAT outputs
  if (is(cat_sim_output, "cat_output")) { # There is a single cat_output element
     nrows <- 1
  } else if (all(sapply(cat_sim_output, is, "cat_output"))) {
    nrows <- length(cat_sim_output)
  } else
    stop("All of the elements of 'cat_sim_output' should be 'cat_output' ",
         "class.")

  if (!is.null(cd) && is.null(cd$ip))
    stop("In order to extract a response data, the cat_design object (cd) should
         have an item pool.")
  if (nrows == 1) {
    # Get estimate history
    eh <- cat_sim_output$est_history
    resp <- sapply(eh, "[[", "resp")
    names(resp) <- sapply(eh, function(x) x$item$id)
  } else {
    if (!is.null(cd)) {
      col_names <- cd$ip$id
    } else {
      eh <- lapply(cat_sim_output, "[[", "est_history") # list of est_history
      col_names <- sort(unique(sapply(do.call("c", eh), function(x) x$item$id)))
    }
    resp <- data.frame(matrix(NA, ncol = length(col_names),
                              nrow = length(cat_sim_output)))
    colnames(resp) <- col_names
    # Create an empty data.frame.
    for (i in 1:length(cat_sim_output)) {
      temp_resp <- get_cat_response_data(cat_sim_output[[i]], cd = cd)
      resp[i, names(temp_resp)] <- temp_resp
    }
  }
  if (remove_na)
    resp <- Filter(f = function(x) !all(is.na(x)), resp)
  if (attach_summary) {
    cs <- summary(cat_sim_output)
    # When there is only one cat output, resp is vector.
    if (nrows == 1) resp <- data.frame(t(resp), check.names = FALSE)
    resp <- cbind(cs, resp)
  }
  return(resp)
}


###############################################################################@
############################# get_cat_administered_items #######################
###############################################################################@
#' Get administered items from a CAT output
#'
#' @description This function returns an item pool object of the
#'   administered items using the items in estimate history. If there is one
#' @param cat_sim_output This is a list object containing elements that are
#' "cat_output" class.
#' @return For \code{cat_output} with only one adaptive test, an
#'   \code{Itempool} class object will be returned. For \code{cat_output} with
#'   more than one adaptive tests, a list of \code{Itempool} class objects will
#'   be returned.
#'
#' @author Emre Gonulates
#'
#' @export
#'
#' @examples
#' cd <- create_cat_design(ip = generate_ip(n = 30), next_item_rule = 'mfi',
#'                         termination_rule = 'max_item',
#'                         termination_par = list(max_item = 10))
#' cat_data <- cat_sim(true_ability = rnorm(10), cd = cd)
#' get_cat_administered_items(cat_data)
get_cat_administered_items <- function(cat_sim_output) {
  if (is(cat_sim_output, "cat_output")) { # There is a single cat_output element
    return(get_administered_items_cpp(cat_sim_output$est_history))
  } else if (all(sapply(cat_sim_output, is, "cat_output"))) {
    eh <- lapply(cat_sim_output, "[[", "est_history") # list of est_history
    return(lapply(eh, function(x) get_administered_items_cpp(x)))
  } else
    stop("All of the elements of 'cat_sim_output' should be 'cat_output' ",
         "class.")
}


###############################################################################@
############################# calculate_exposure_rates #########################
###############################################################################@
#' Calculate exposure rate of items for CAT
#' @description This function calculates the exposure rate of items for a
#' CAT. It takes a list of \code{cat_output} objects and \code{cat_design}
#' object and returns exposure rate of each item.
#'
#' @param cat_sim_output This is a list object containing elements that are
#' "cat_output" class.
#' @param cd A \code{cat_design} object that is created by function
#'          \code{create_cat_design}.
#' @param item_ids A vector of Item (or Testlet) ids in the item pool.
#' @return This function returns a numeric vector of each item's exposure rate
#'   where the names of each exposure rate value is the item's id.
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @seealso \code{\link{cat_sim}}
#'
#' @examples
#' cd <- create_cat_design(ip = generate_ip(n = 30), next_item_rule = 'mfi',
#'                         termination_rule = 'max_item',
#'                         termination_par = list(max_item = 10))
#' cat_data <- cat_sim(true_ability = rnorm(10), cd = cd)
#' calculate_exposure_rates(cat_data, cd = cd)
#'
calculate_exposure_rates <- function(cat_sim_output, cd = NULL, item_ids = NULL) {
  # Check whether it is one or multiple CAT outputs
  if (is(cat_sim_output, "cat_output")) { # There is a single cat_output element
    cat_sim_output <- list(cat_sim_output)
  # if it is not a list of cat_output
  } else if (!is.list(cat_sim_output) ||
             !all(sapply(cat_sim_output, is, "cat_output")))
    stop("All of the elements of 'cat_output' should be 'cat_output' class.")

  if (is.null(item_ids)) {
    if (is.null(cd)) {
      stop("Either 'cd' or 'item_ids' should be provided.")
    } else {
      if (is.null(cd$ip))
        stop("In order to calculate exposure rates, the cat_design object (cd)
              should have an item pool.")
      item_ids <- get_slot_itempool_cpp(cd$ip, "id")
    }
  }

  return(calculate_exposure_rates_cpp(item_ids, cat_sim_output))
}


###############################################################################@
############################# calculate_overlap_rates ##########################
###############################################################################@
#' Calculate overlap rate of items for CAT
#' @description This function calculates the overlap rate of items for a
#' CAT. It takes a list of \code{cat_output} objects and \code{cat_design}
#' object and returns exposure rate of each item.
#'
#' @param cat_sim_output This is a list object containing elements that are
#' "cat_output" class.
#' @param cd A \code{cat_design} object that is created by function
#'          \code{create_cat_design}.
#' @param item_ids A vector of item (or Testlet) ids in the item pool.
#' @return This function returns a numeric vector of each item's overlap rate
#'   where the names of each overlap rate value is the item's id.
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @seealso \code{\link{cat_sim}}
#'
#' @examples
#' cd <- create_cat_design(ip = generate_ip(n = 30), next_item_rule = 'mfi',
#'                         termination_rule = 'max_item',
#'                         termination_par = list(max_item = 10))
#' cat_data <- cat_sim(true_ability = rnorm(10), cd = cd)
#' calculate_overlap_rates(cat_data, cd = cd)
#'
calculate_overlap_rates <- function(cat_sim_output, cd = NULL, item_ids = NULL) {
  # Check whether it is one or multiple CAT outputs
  if (is(cat_sim_output, "cat_output")) { # There is a single cat_output element
    cat_sim_output <- list(cat_sim_output)
  # if it is not a list of cat_output
  } else if (!is.list(cat_sim_output) ||
             !all(sapply(cat_sim_output, is, "cat_output")))
    stop("All of the elements of 'cat_output' should be 'cat_output' class.")
  if (is.null(item_ids)) {
    if (is.null(cd)) {
      stop("Either 'cd' or 'item_ids' should be provided.")
    } else {
      if (is.null(cd$ip))
        stop("In order to calculate overlap rates, the cat_design object (cd)
              should have an item pool.")
      item_ids <- get_slot_itempool_cpp(cd$ip, "id")
    }
  }
  return(calculate_overlap_rates_cpp(item_ids, cat_sim_output))
}
