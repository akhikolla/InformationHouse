
#' Item Parameter Drift
#'
#' @description This function detects the unstable (i.e. items whose item
#'   parameter values drifted) for a given two sets of items.
#'
#' @param ip1 Itempool object for the first calibration.
#' @param ip2 Itempool object for the second calibration.
#' @param method The method of item parameter drift analysis.
#' @param anchor_item_ids Anchor item ids. If \code{NULL}, all items are
#'   assumed to be anchor items.
#' @param alpha Two tailed critical value to detect the unstable items. For
#'   example if \deqn{alpha = 0.05}, the critical value is calculated using
#'   \code{qnorm(1-alpha/2)} (= 1.96). Items whose absolute robust z values
#'   are larger than this number will be flagged as unstable.
#'
#'   \describe{
#'     \item{"robust-z"}{Robust-Z method based on the Huynh and Meyer (2010). }
#'   }
#'
#' @return Return a list depending on the method:
#'   \describe{
#'     \item{robust-z}{
#'       \describe{
#'         \item{\code{output$a$cor}}{Correlation between two $a$ parameter
#'           sets.}
#'         \item{\code{output$a$sd_ratio}}{The ratio of the standard deviation
#'           of \code{ip2} to the standard deviation of \code{ip1}.}
#'         \item{\code{output$a$robust_z}}{Robust-z statistic values for each
#'           item discrimination parameter.}
#'         \item{\code{output$a$unstable}}{Item id's which were flagged if
#'           robust z statistic value for a parameters is larger than the
#'           absolute value of  the critical value
#'           (i.e. \code{qnorm(1-alphe/2)}).}
#'         \item{\code{output$b$robust_z}}{Robust-z statistic values for each
#'           item difficulty or threshold parameter. If an item has threshold
#'           parameters, robust z statistic will be calculated for each
#'           threshold.}
#'         \item{\code{output$b$unstable}}{Item id's which were flagged if
#'           robust z statistic for difficulty/threshold parameters are larger
#'           than the absolute value of the critical value (i.e.
#'           \code{qnorm(1-alphe/2)}).}
#'       }
#'     }
#'   }
#'
#'
#' @author Emre Gonulates
#'
#' @export
#'
#' @references
#' Huynh, Huynh and Meyer, Patrick (2010) "Use of Robust z in Detecting
#' Unstable Items in Item Response Theory Models,"
#' \emph{Practical Assessment, Research, and Evaluation}: Vol. 15 , Article 2.
#' DOI: \url{https://doi.org/10.7275/ycx6-e864}
#' Available at: \url{https://scholarworks.umass.edu/pare/vol15/iss1/2/}
#'
#' @examples
#' # The example from Huynh and Meyer (2010)
#' ip1 <- c(itempool(
#'   a = c(0.729, 0.846, 0.909, 0.818, 0.742, 0.890, 1.741, 0.907, 1.487, 1.228,
#'         0.672, 1.007, 1.016, 0.776, 0.921, 0.550, 0.624, 0.984, 0.506, 0.594,
#'         0.687, 0.541, 0.691, 0.843, 0.530, 0.462, 1.007, 0.825, 0.608, 1.177,
#'         0.900, 0.861, 0.843, 1.404, 0.446, 1.014, 1.632, 0.831, 1.560, 0.798),
#'   b = c(1.585, 0.635, -0.378, -0.100, -0.195, 0.749, 1.246, 1.016, -0.234,
#'         0.537, 0.070, 1.985, 1.101, -0.742, 0.463, -0.060, 0.477, 1.084,
#'         -2.340, 1.068, -0.055, -1.045, 1.859, 0.645, -0.689, -2.583, 1.922,
#'         0.709, 0.499, 1.973, 0.104, 0.809, 0.640, 0.247, 0.820, 1.837,
#'         2.129, 1.012, 1.774, 0.095),
#'   c = c(0.134, 0.304, 0.267, 0.176, 0.215, 0.194, 0.267, 0.159, 0.095,
#'         0.197, 0.089, 0.272, 0.229, 0.159, 0.162, 0.100, 0.259, 0.167,
#'         0.000, 0.242, 0.323, 0.000, 0.196, 0.189, 0.000, 0.000, 0.334,
#'         0.538, 0.125, 0.511, 0.192, 0.353, 0.103, 0.241, 0.245, 0.118,
#'         0.155, 0.132, 0.215, 0.148),
#'   model = "3PL"),
#'   item(a = 0.561, b = c(0.784, -0.113, 1.166), model = "GPCM"),
#'   item(a = 0.745, b = c(3.687, 2.506, -0.001), model = "GPCM"))
#'
#' ip2 <- c(itempool(
#'   a = c(0.650, 0.782, 0.816, 0.787, 0.611, 0.888, 1.192, 0.589, 1.211,
#'         0.742, 0.526, 0.690, 0.996, 0.816, 0.781, 0.507, 0.378, 0.976,
#'         0.473, 0.364, 0.585, 0.566, 0.511, 0.718, 0.354, 1.080, 0.840,
#'         0.865, 0.528, 0.814, 0.555, 0.701, 0.530, 1.220, 0.344, 0.966,
#'         1.044, 0.358, 1.192, 0.615),
#'   b = c(0.676, -0.525, -1.749, -1.092, -1.619, -0.406, -0.132, 0.006,
#'         -1.352, -0.872, -1.242, 0.873, 0.239, -2.038, -0.487, -1.372,
#'         -1.492, 0.214, -4.537, 0.220, -0.686, -2.394, 0.747, -0.467,
#'         -3.629, -5.000, 0.927, 0.305, -0.839, 1.270, -1.618, -0.091,
#'         -1.228, -1.019, -1.453, 1.090, 1.743, -1.436, 1.024, -1.358),
#'   c = c(0.110, 0.316, 0.161, 0.149, 0.145, 0.200, 0.243, 0.059, 0.081,
#'         0.075, 0.028, 0.267, 0.242, 0.189, 0.184, 0.121, 0.000, 0.170,
#'         0.000, 0.151, 0.383, 0.000, 0.195, 0.177, 0.000, 0.000, 0.352,
#'         0.647, 0.116, 0.501, 0.000, 0.286, 0.000, 0.248, 0.064, 0.150,
#'         0.126, 0.000, 0.187, 0.007),
#'   model = "3PL"),
#'   item(a = 0.486, b = c(-0.539, -1.489, -0.052), model = "GPCM"),
#'   item(a = 0.737, b = c(2.599, 1.250, -1.209), model = "GPCM"))
#' ipd(ip1, ip2)
#'
ipd <- function(ip1, ip2, method = "robust-z", anchor_item_ids = NULL,
                alpha = 0.01) {
  if (method == "robust-z") {
    # Check if anchor item ids is not null. If it is not, pull those items
    if (!is.null(anchor_item_ids)) {
      # Check if anchor_item_ids appeared in both item pools. If not raise an
      # error.
      if (!all(anchor_item_ids %in% ip2$resp_id) ||
          !all(anchor_item_ids %in% ip1$resp_id))
        stop("All 'anchor_item_ids' should be in both ip2 and ip1.")
      ip2 <- ip2[anchor_item_ids %in% ip2$resp_id]
      ip1 <- ip1[anchor_item_ids %in% ip1$resp_id]
    }

    output <- list()
    a2 <- ip2$a
    a1 <- ip1$a
    if ((is.null(a2) || any(is.na(a2))) ||
        (is.null(a1) || any(is.na(a1))))
      stop("In the item pools, there are some items that do not have 'a' ",
           "parameters.")
    a_diff <- log(a2) - log(a1)
    output$a$robust_z <- (a_diff - stats::median(a_diff)) /
      (0.74 * stats::IQR(a_diff))
    output$a$cor <- stats::cor(a2, a1)
    output$a$sd_ratio <- stats::sd(a2)/stats::sd(a1)

    cv <- stats::qnorm(1 - alpha/2) # two tailed critical value
    output$a$unstable <- ip1$id[abs(output$a$robust_z) >= cv]

    # Select only stable items to calculate the linking constant
    A <- exp(mean(a_diff[abs(output$a$robust_z) < cv]))

    ip1_list <- flatten_itempool_cpp(ip1)
    ip2_list <- flatten_itempool_cpp(ip2)

    b1 <- b2 <- c()
    for (i in seq_along(ip1_list)) {
      temp_b <- ip1_list[[i]]@parameters$b
      b1 <- c(b1, stats::setNames(
        temp_b, paste0(ip1_list[[i]]$id, if (length(temp_b) > 1)
          paste0(".", 1:length(temp_b)) else NULL)))

      temp_b <- ip2_list[[i]]@parameters$b
      b2 <- c(b2, stats::setNames(
        temp_b, paste0(ip1_list[[i]]$id, if (length(temp_b) > 1)
          paste0(".", 1:length(temp_b)) else NULL)))

    }
    b_diff <- b1 - A * b2
    output$b$robust_z <- (b_diff - stats::median(b_diff)) /
      (0.74 * stats::IQR(b_diff))
    output$b$unstable <- names(output$b$robust_z)[output$b$robust_z >= cv]
    return(output)
  } else
    stop("This method has not been implemented yet.", call. = FALSE)
}
