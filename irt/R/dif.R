
###############################################################################@
############################# dif ##############################################
###############################################################################@
#' Evaluate Differential Item Functioning (DIF) of a test
#'
#' @description
#' \code{dif} evaluates Differential Item Functioning (DIF) of a test.
#'
#'
#' @param resp A vector of item responses.
#' @param group Group membership
#' @param focal_name In the group variable, the value that represents the focal
#'   group.
#' @param ip An \code{\link{Itempool-class}} object.
#' @param type The type of the DIF method.
#'
#' @return A data.frame of DIF values.
#'
#' @include Item-class.R
#' @include Itempool-class.R
#' @include Item-class-methods.R
#' @include Itempool-class-methods.R
#'
#' @author Emre Gonulates
#'
#' @export
#'
dif <- function(resp, group, focal_name, ip = NULL, type = "mh") {
  if (type == "mh") {
    result <- data.frame(id = paste0("Item", 1:ncol(resp)), mh_statistic = NA,
                         chisq = NA, p_value = NA, ETS = NA, ETS_class= NA)
    if (!is.null(colnames(resp))) result$id <- colnames(resp)
    total_score <- rowSums(resp, na.rm = TRUE)
    K <- sort(unique(total_score))
    group[group == focal_name] <- "focal"
    group[group != "focal"] <- "reference"
    group <- as.factor(group)

    for (i in seq_len(ncol(resp))) { # i = item
      n_r1 <- rep(NA, length(K))
      n_r0 <- rep(NA, length(K))
      n_f1 <- rep(NA, length(K))
      n_f0 <- rep(NA, length(K))
      n_r <- rep(NA, length(K))
      n_f <- rep(NA, length(K))
      n_0 <- rep(NA, length(K))
      n_1 <- rep(NA, length(K))
      n <- rep(NA, length(K))
      for (k in seq_len(length(K))) {
        temp <- which(total_score == K[k])
        temp_resp <- factor(resp[temp, i], levels = unique(resp[, i]))
        temp <- stats::addmargins(table(group[temp], temp_resp))
        n_r1[k] <- temp["reference", "1"]
        n_r0[k] <- temp["reference", "0"]
        n_f1[k] <- temp["focal", "1"]
        n_f0[k] <- temp["focal", "0"]
        n_r[k] <- temp["reference", "Sum"]
        n_f[k] <- temp["focal", "Sum"]
        n_0[k] <- temp["Sum", "0"]
        n_1[k] <- temp["Sum", "1"]
        n[k] <- temp["Sum", "Sum"]
      }
      temp <- n > 1
      result$mh_statistic[i] <- sum((n_r1 * n_f0 / n)[temp]) / sum((n_r0 * n_f1 / n)[temp])
      mu <- (n_r * n_1 / n)[temp]
      v <- (n_r * n_f * n_1 * n_0 / (n^2*(n-1)))[temp]
      result$chisq[i] <- (abs(sum(n_r1[temp]) - sum(mu)) - 0.5)^2 / sum(v)
    }
    result$ETS <- -2.35 * log(result$mh_statistic)
    result$p_value <- 1 - stats::pchisq(q = result$chisq, df = 1)
    result$ETS_class <- ifelse(abs(result$ETS) < 1, "A",
                               ifelse(abs(result$ETS) < 1.5, "B", "C"))
  }
  return(result)
}
