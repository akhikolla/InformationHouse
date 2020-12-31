
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% rsss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Convert raw score to scale score and vice versa
#'
#' @param ip An \code{\link{Itempool-class}} object.
#' @param raw_score A value (or vector of values) representing raw score(s).
#' @param scale_score A value (or vector of values) representing scale score(s).
#' @param theta_range The limits of the scale score. The default is
#'   \code{c(-5, 5)}.
#'
#' @return A vector of raw or scale scores.
#'
#'
#' @importFrom stats uniroot
#' @export
#'
#' @author Emre Gonulates
#'
rsss <- function(ip, raw_score = NULL, scale_score = NULL,
                 theta_range = c(-5, 5)) {
  result <- NULL
  if (!is.null(raw_score)) {
    max_possible_score <- sum(ip$resp_max_score)
    if (any(raw_score > max_possible_score) || any(raw_score < 0))
      stop(paste0("Raw score values cannot larger than the maximum possible ",
                  "raw score of ", max_possible_score, " or smaller than 0."))
    raw_score_to_scale_score <- Vectorize(function(y) {
      uniroot(f = function(x) {
        sum(prob(ip = ip, theta = x, expected_value = TRUE)) - y
        }, lower=sort(theta_range)[1], upper = sort(theta_range)[2])$root
    })
    result <- tryCatch(
      raw_score_to_scale_score(raw_score),
      error=function(e) {
        if (grepl(pattern = "values at end points not of opposite sign",
                  e$message))
          stop(paste0(
            "\nScale score cannot be calculated for some raw scores. Please ",
            "provide a wider 'theta_range' than c(", theta_range[1], " ,",
            theta_range[2], "). \nOr, remove raw score(s): ",
            paste0(c(switch(0 %in% raw_score, 0, NULL),
                     switch(max_possible_score %in% raw_score,
                            max_possible_score, NULL)), collapse = ", "), "."),
               call. = FALSE)
      }
    )
  } else if (!is.null(scale_score) && all(sapply(scale_score, is.numeric))) {
    scale_score_to_raw_score <- Vectorize(function(x) {
      sum(prob(ip = ip, theta = x, expected_value = TRUE))
    })
    result <- scale_score_to_raw_score(x = scale_score)
  }
  return(result)
}

