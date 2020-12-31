# Core functionality ------------------------------------------------------

#' An S4 class containing summary information about a \code{transition} object
#'
#' @slot result a data frame with the summary information
summaryTransition <- setClass(
  "summaryTransition", representation(
    result = "data.frame", trans_object = "list"
  )
)

#' Evaluate the effectiveness of predicted transitions for an object of class
#' \code{transition}
#'
#' @param object a \code{transition} object to analyze
#' @inheritParams get_transition_info
#' @param ... further arguments passed to or from methods, currently unused
#'
#' @return a data frame containing indicators that reflect, in different ways,
#'   the effectiveness of predicted transitions compared to the set of actual
#'   (reference) transitions
#' @export
#'
#' @examples
#' predictions <- (sample(1:100)%%2)
#' references  <- (sample(1:100)%%2)
#' window_size <- 7
#' transitions <- get_transition_info(predictions, references, window_size)
#' summary(transitions)
summary.transition <- function(
  object, ...
) {

  # Marginal totals

    reference_positives <- sum(object$references)
    predicted_positives <- sum(object$predictions)

  # Cell totals

    rejects <- object$matchings$rejected
    true_positives <- sum(!rejects)

  # Lags & RMSE

    abs_lags <-
      object$matchings$abs_lag[!rejects] %>%
      mean_sd(digits = 1, nsmall = 1) %>%
      sapply(
        function(x) if(is.nan(x)) NA_real_ else x,
        simplify = FALSE
      ) %>%
      c(stringsAsFactors = FALSE, row.names = list(NULL)) %>%
      do.call(data.frame, .)

    signed_lags <-
      object$matchings$signed_lag[!rejects] %>%
      mean_sd(digits = 1, nsmall = 1) %>%
      sapply(
        function(x) if(is.nan(x)) NA_real_ else x,
        simplify = FALSE
      ) %>%
      c(stringsAsFactors = FALSE, row.names = list(NULL)) %>%
      do.call(data.frame, .)

    rmse <-
      object$matchings$abs_lag[!rejects]^2 %>%
      mean(.) %>%
      sqrt(.) %>%
      round(1) %>%
      sapply(function(x) if(is.nan(x)) NA_real_ else x)

    rmse_prop <- 1 - rmse/object$window_size

  # Recall/precision

    recall <- true_positives / reference_positives

    precision <- true_positives / predicted_positives

  # Summarize

    class(object) <- "list"
    data.frame(

      window_size = object$window_size,

      reference_positives = reference_positives,
      predicted_positives = predicted_positives,

      true_positives = true_positives,

      n_rejected_pairs = sum(rejects),

      mean_abs_lag_indices = abs_lags$mean,
      sd_abs_lag_indices = abs_lags$sd,
      mean_sd_abs_lag_indices = gsub(
        "NaN", "NA", abs_lags$sum_string
      ),

      mean_signed_lag_indices = signed_lags$mean,
      sd_signed_lag_indices = signed_lags$sd,
      mean_sd_signed_lag_indices = gsub(
        "NaN", "NA", signed_lags$sum_string
      ),

      recall = recall,
      precision = precision,

      rmse_lag_indices = rmse,
      rmse_prop = rmse_prop,

      row.names = NULL,
      stringsAsFactors = FALSE

    ) %>% cbind(
      ., aggregated_performance = rowMeans(
        .[ ,c("recall", "precision", "rmse_prop")]
      )
    ) %>%
    new("summaryTransition", result = ., trans_object = object)

}

# Addition and Subtraction ------------------------------------------------

#' Addition and subtraction for objects of class \code{summaryTransition}
#'
#' @param e1 the first object
#' @param e2 the second object
#'
#' @keywords internal
add_summaryTransition <- function(e1, e2) {

  not_used <- "Not used in combined objects"

  e1 <- e1@trans_object
  e2 <- e2@trans_object

  references <- c(
    e1$references, e2$references
  )

  predictions <- c(
    e1$predictions, e2$predictions
  )

  matchings <- rbind(
    e1$matchings, e2$matchings
  )

  window_size <-
    e1$window_size %>%
    c(e2$window_size) %>%
    unique() %T>%
    {if (length(.) > 1) stop(
      "Cannot add if spurious pairing thresholds",
      " differ for the objects",
      call. = FALSE
    )} %>%
    as.numeric()

  list(
    window_size = window_size,
    references = references,
    predictions = predictions,
    reference_transition_indices = not_used,
    prediction_transition_indices = not_used,
    pruned_reference_transition_indices = not_used,
    pruned_prediction_transition_indices = not_used,
    false_negative_indices = not_used,
    false_positive_indices = not_used,
    reference_preferences = not_used,
    prediction_preferences = not_used,
    matchings = matchings
  ) %>%
  summary.transition(.)

}

#' @rdname add_summaryTransition
#' @export
setMethod("+", signature(e1 = "summaryTransition", e2 = "summaryTransition"),
  function(e1, e2) add_summaryTransition(e1, e2))

#' @rdname add_summaryTransition
#' @keywords internal
subtract_summaryTransition <- function(e1, e2) {

  diff_vars <- c(
    "window_size","recall", "precision",
    "mean_abs_lag_indices", "mean_signed_lag_indices",
    "rmse_lag_indices", "rmse_prop", "aggregated_performance"
  )

  result <-
    sapply(
      diff_vars,
      function(x) e1@result[ ,x] - e2@result[ ,x],
      simplify = FALSE
    ) %>%
    do.call(data.frame, .) %>%
    {stats::setNames(., paste("diff", names(.), sep = "_"))} %>%
    data.frame(
      window_size = e1@result$window_size, .,
      stringsAsFactors = FALSE, row.names = NULL
    )

  if(result$diff_window_size != 0) {
    warning(
      "Spurious pairing thresholds differ for the objects",
      call. = FALSE
    )
    result$window_size <- paste(
      e1@result$window_size, e2@result$window_size, sep = "--"
    )
  }

  list(differences = result, e1 = e1@result, e2 = e2@result)

}

#' @rdname add_summaryTransition
#' @export
setMethod("-", signature(e1 = "summaryTransition", e2 = "summaryTransition"),
  function(e1, e2) subtract_summaryTransition(e1, e2))

# Coercion ----------------------------------------------------------------

#' As("summaryTransition", "data.frame")
#'
#' @name as
#' @family summaryTransition
setAs("summaryTransition", "data.frame", function(from) from@result)

#' As("summaryTransition", "list")
#'
#' @name as
#' @family summaryTransition
setAs("summaryTransition", "list", function(from) as.list(from@result))
