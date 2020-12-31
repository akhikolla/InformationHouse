

#' Item Analysis Function
#'
#' @param resp A \code{matrix} or \code{data.frame} containing the item
#'   responses.
#' @param criterion Provide a continuous criterion variable such as a total
#'   raw score, or theta score that will be used in the calculation of
#'   correlation calculations. If this value is \code{NULL}, the total score
#'   will be used.
#' @param suppress_output If \code{TRUE}, the function will suppress
#'   console output. Default value is \code{FALSE}
#'
#' @return A list of
#'   \describe{
#'     \item{'id'}{Item ID.}
#'     \item{'n'}{Number of examinees responded this item.}
#'     \item{'pval'}{p-value, proportion of examinees correctly answered items.}
#'     \item{'pbis'}{Point biserial correlation.}
#'     \item{'bis'}{Biserial correlation.}
#'     \item{'pbis_adj'}{Point biserial correlation between item and total score
#'       without this item.}
#'     \item{'bis_adj'}{Biserial correlation between item and total score
#'       without this item.}
#'   }
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @examples
#' theta <- rnorm(100)
#' ip <- generate_ip(n = 20)
#' resp <- sim_resp(ip = ip, theta = theta, prop_missing = .2)
#' # Item analysis based on total scores
#' item_analysis(resp)
#' # Item analysis based on theta scores
#' item_analysis(resp, criterion = theta)
#'
item_analysis <- function(resp, criterion = NULL, suppress_output = FALSE) {
  # Calculate the number of non-missing responses
  output <- data.frame(row.names = paste0(1:ncol(resp)))
  # Add item ids
  if (is.null(colnames(resp))) output$id <- rep(NA, ncol(resp)) else
    output$id <- colnames(resp)
  # Add number of
  output$n <- colSums(!is.na(resp))
  # Calculate the p-value
  output$pval <- colMeans(resp, na.rm = TRUE)

  if (is.null(criterion)) theta <- rowSums(resp, na.rm = TRUE) else
    theta <- criterion
  # Calculate point-biserial correlation
  output$pbis <- apply(resp, 2, biserial, theta, "point-biserial")
  # Calculate biserial correlation
  output$bis <- apply(resp, 2, biserial, theta, "default")

  if (is.null(criterion)) {
    # Calculate corrected point-biserial correlation, i.e. remove the item from
    # the calculation of the total score.
    output$pbis_adj <- sapply(1:ncol(resp), function(i)
      biserial(score = resp[, i], total_score = rowSums(resp[, -i], na.rm = TRUE),
               method = "point-biserial"))
    output$bis_adj <- sapply(1:ncol(resp), function(i)
      biserial(score = resp[, i], total_score = rowSums(resp[, -i], na.rm = TRUE),
               method = "default"))
  }
  return(output)
}


#' Calculate point-biserial correlation
#'
#' @param score Item scores of each examinee for which point-biserial
#'   correlation will be calculated
#' @param total_score Total score of each examinee
#'
#' @return Point-biserial correlation value
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
point_biserial <- function(score, total_score) {
  # Only use the complete observations
  complete_obs <- !(is.na(score) | is.na(total_score))
  score <- score[complete_obs]
  total_score <- total_score[complete_obs]
  # Mean of total scores of  examinees who correctly answered the item
  mu_correct <- mean(total_score[score == 1], na.rm = TRUE)
  mu_all <- mean(total_score, na.rm = TRUE)
  # Item difficulty
  p <- mean(score, na.rm = TRUE)
  # Use n instead of (n-1) for the calculation of standard deviation
  # sigma <- stats::sd(total_score, na.rm = TRUE)
  sigma <- sqrt(sum((total_score - mean(total_score))^2)/length(total_score))
  return((mu_correct - mu_all) * sqrt(p/(1-p)) / sigma)
}

#' Calculate biserial correlation
#'
#' @param score Item scores of each examinee for which biserial correlation will
#'   be calculated
#' @param total_score Total score of each examinee
#' @param method Type of the biserial correlation calculation method.
#'   \describe{
#'     \item{\strong{"default"}}{The most common way to calculate biserial
#'     correlation. }
#'     \item{\strong{"point-biserial"}}{Calculate point-biserial correlation. }
#'     \item{\strong{"clemans-lord"}}{Modified biserial correlation value based
#'       on Clemans (1958) and Lord (1962).}
#'     \item{\strong{"brogden"}}{Modified biserial correlation value based on:
#'       Brogden (1949)}
#'     \item{\strong{"rank"}}{Rank biserial correlation value based on Cureton
#'       (1968).}
#'     }
#'
#' @return Biserial correlation value
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @references
#' Brogden, H. E. (1949). A new coefficient: Application to biserial
#'   correlation and to estimation of selective efficiency. Psychometrika, 14,
#'   169-182.
#'
#' Clemans, W. V. (1958) An index of item-criterion relationship. Educational
#'   and Psychological Measurement, 18, 167-172.
#'
#' Cureton, E. E. (1968). Rank biserial correlation when ties are present.
#'   Educational and Psychological Measurement, 28, 77-79.
#'
#' Kraemer, H. C. (1981). Modified biserial correlation coefficients.
#'   Psychometrika, 46(3), 275-282.
#'
#' Lord, F. M. (1963). Biserial estimates of correlation. Psychometrika, 28,
#'   81–85.
#'
#' @examples
#' # The example is from Salkind, Rasmussen (2007) Encyclopedia of measurement
#' # and statistics, pages 94-97
#' score <- c(rep(0, 16), rep(1, 22))
#' total_score <- c(87, 90, 94, 94, 97, 103, 103, 104, 106, 108, 109, 109, 109,
#'                  112, 119, 132, 100, 103, 103, 106, 112, 113, 114, 114, 118,
#'                  119, 120, 120, 124, 133, 135, 135, 136, 141, 155, 157, 159,
#'                  162)
#' # Calculate biserial correlation
#' biserial(score, total_score)
#' # Calculate point-biserial correlation
#' biserial(score, total_score, method = "point-biserial")
#' # Calculate modified biserial correlation (based on Brogden (1949))
#' biserial(score, total_score, method = "brogden")
#' # Calculate modified biserial correlation (Clemans-Lord)
#' biserial(score, total_score, method = "clemans-lord")
#'
biserial <- function(score, total_score, method = "default") {
  if (method == "point-biserial")
    return(stats::cor(score, total_score, use = "pairwise.complete.obs"))

  # Only use the complete observations
  complete_obs <- !(is.na(score) | is.na(total_score))
  score <- score[complete_obs]
  total_score <- total_score[complete_obs]

  # n <- sum(!is.na(score))
  n <- as.double(length(score))

  if (method == "rank") {
    # Kraemer pointed out that Cureton's (1958, 1964) rank biserial correlation
    # coefficient "essentially replaces observations with their ranks and then
    # applies Brogden's approach. " (p.280)
    total_score <- rank(total_score)
    method <- "brogden"
  }
  if (method == "clemans-lord") {
    dev <- total_score - mean(total_score)
    num <- sum(score * dev)
    if (num < 0)
      return(-(num / sum((1 - sort(score, decreasing = TRUE)) * dev)))
    # For positive values Lord's modification is the same as Brogden.
    return(sum(score * dev) /  sum(sort(score, decreasing = TRUE) * dev))
  }
  if (method == "brogden") {
    dev <- total_score - mean(total_score)
    return(sum(score * dev) /  sum(sort(score, decreasing = TRUE) * dev))
  }


  # Mean of total scores of  examinees who correctly answered the item
  mu1 <- mean(total_score[score == 1], na.rm = TRUE)
  mu0 <- mean(total_score[score == 0], na.rm = TRUE)
  n1 <- as.double(sum(score == 1, na.rm = TRUE))
  n0 <- as.double(sum(score == 0, na.rm = TRUE))
  # # The following is also valid calculation but the above one is shorter.
  # if (method == "brogden") {
  #   ranked <- sort(total_score, decreasing = TRUE)
  #   D <-  mean(ranked[1:n1]) - mean(ranked[(n1+1):(n0+n1)])
  #   return((mu1 - mu0) / D)
  # }
  sigma <- stats::sd(total_score, na.rm = TRUE)
  u <- stats::dnorm(stats::qnorm(n1/n))
  return(((mu1 - mu0) / sigma) * (n1 * n0 / (u * n^2)))
}



#' Distractor Analysis Function
#'
#' @param raw_resp A \code{matrix} or \code{data.frame} containing the item
#'   responses.
#' @param key The answer key for the responses
#' @param criterion Provide a continuous criterion variable such as a total
#'   raw score, or theta score that will be used in the calculation of
#'   correlation calculations. If this value is \code{NULL}, the total score
#'   will be used.
#' @param adjusted If \code{TRUE}, the biserial will be calculated using the
#'   total score that is calculated using all of the items except the item
#'   that biserial is calculated for.
#' @param suppress_output If \code{TRUE}, the function will suppress
#'   console output. Default value is \code{FALSE}
#'
#' @return A list of
#'   \describe{
#'     \item{'prop'}{Observed proportions of each choice.}
#'     \item{'biserial'}{Biserial correlation between the examinees selected
#'           the choice and the total scores.}
#'     \item{'score'}{Scored response matrix.}
#'   }
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @examples
#' n_item <- 10 # sample(8:12, 1)
#' n_theta <- 50 # sample(100:200, 1)
#' raw_resp <- matrix(sample(LETTERS[1:4], n_item * n_theta, replace = TRUE),
#'                    nrow = n_theta, ncol = n_item,
#'                    dimnames = list(paste0("Examinee-", 1:n_theta),
#'                                    paste0("Item-", 1:n_item)))
#' # Add some missing responses
#' raw_resp[sample(1:length(raw_resp), round(length(raw_resp)*.1))] <- NA
#' # Prepare answer key
#' key <- sample(LETTERS[1:4], n_item, replace = TRUE)
#'
#' # Run distractor analysis:
#' da <- distractor_analysis(raw_resp = raw_resp, key = key)
#'
distractor_analysis <- function(raw_resp, key, criterion = NULL,
                                adjusted = FALSE, suppress_output = FALSE) {
  # possible_options <- unique(as.vector(unlist(raw_resp))) # Very very slow
  # possible_options <- unique(as.vector(apply(raw_resp, 2, unique))) # Slow
  if (inherits(raw_resp, "matrix")) raw_resp <- as.data.frame(raw_resp)
  #   stop("Invalid 'raw_resp'. Please provide a matrix or a data.frame")
  possible_options <- unique(rapply(raw_resp, unique, how = "unlist"))
  possible_options <- sort(possible_options[!is.na(possible_options)])

  ### Calculate the observed proportion of each choice option
  choice_props <- t(sapply(lapply(raw_resp, factor, levels = possible_options),
                           function(x) prop.table(table(x))))
  choice_props <- cbind.data.frame(n = colSums(!is.na(raw_resp)), choice_props)

  # ### Calculate the observed proportion of each choice option
  # choice_props <- t(apply(raw_resp, 2, function(x) prop.table(table(x))))
  # # Add total number or examinees answered that item:
  # choice_props <- cbind.data.frame(
  #   n = apply(raw_resp, 2, function(x) sum(!is.na(x))), choice_props)

  ### Calculate the biserial correlation of each choice option.
  # Calculate the score matrix
  resp <- score_raw_resp(raw_resp, key)
  choice_bis <- data.frame(key = key, row.names = colnames(raw_resp))

  # If not adjusted use the same total_score for all items
  if (is.null(criterion)) total_score <- rowSums(resp, na.rm = TRUE) else
    total_score <- criterion
  if (adjusted && is.null(criterion)) {
    total_score <- matrix(total_score, ncol = ncol(resp), nrow = nrow(resp))
    total_score <- ifelse(is.na(raw_resp), total_score, total_score - resp)
    for (po in possible_options)
      choice_bis[[po]] <- sapply(1:ncol(raw_resp), function(i) biserial(
        score = (raw_resp[, i] == po) * 1, total_score = total_score[, i]))
  } else {
    for (po in possible_options)
      choice_bis[[po]] <- apply((raw_resp == po) * 1, 2, biserial, total_score)
  }
  output <- list(prop = choice_props, biserial = choice_bis)
  if (!suppress_output) print(output)
  output$scores <- resp

  return(invisible(output))
}



#' Score Raw Responses
#'
#' @param raw_resp A \code{matrix} or \code{data.frame} containing the item
#'   responses.
#' @param key The answer key for the responses
#'
#' @return Scored response matrix
#'
#' @keywords internal
#'
#' @author Emre Gonulates
#'
score_raw_resp <- function(raw_resp, key) {
  output <- t(apply(raw_resp, 1, function(x) as.integer(x == key)))
  colnames(output) <- colnames(raw_resp)
  return(output)
}
