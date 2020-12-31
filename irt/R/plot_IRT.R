
############################################################################@###
############################# plot.Item ####################################@###
############################################################################@###
#' Plot Item Characteristic Curve of an \code{Item} object
#' @description
#' \code{plot.Item} Plots the item characteristic curve.
#'
#' @param x An \code{\link{Item-class}} object.
#' @param theta_range The boundaries of x axis.
#' @param title Title of the plot. By default if the item is 1-4PM IRT model
#'          then the title will be "Item Characteristic Curve" if the item
#'          follows Graded Response Model the title will be
#'          "Category Response Functions". Set it \code{NULL} to remove it.
#' @param suppress_plot If \code{FALSE} the function will print the plot. If
#'          \code{TRUE}, function will return the plot object. Default value is
#'          \code{FALSE}.
#' @param category_names If the model used is 'GRM' (Graded Response Model)
#'          these names will serve as category names. For example,
#'          c("Strongly Disagree", "Disagree", "Agree", "Strongly Agree").
#'          The default is \code{FALSE} where the default category scores
#'          will be printed. If the value is \code{NULL} no legend will be
#'          printed but the categories will be printed differently.
#' @param legend_title The title of the plot's legend.
#' @param ... Additional arguments that will be passed to \code{geom_line}
#'
#' @return Depending on the value of \code{suppress_plot} function either prints
#' the item characteristic curve or returns the plot object.
#'
#' @export
#' @importFrom ggplot2 aes aes_string element_text geom_line ggplot ggtitle
#'             guide_legend guides theme theme_bw xlab ylab ylim
#'             scale_x_continuous scale_color_discrete
#' @importFrom stats reshape runif
#'
#' @author Emre Gonulates
#'
#' @examples
#' plot(x = item(b = 0.3, D = 1))
#'
#' item <- item(a = 1.2, b = 0.3, c = .2)
#' plot(item)
#' plot(item(a = 1.2, b = 0.3, c = .2, d = .89, D = 1))
#'
#' # Plot Graded Response Model
#' ip <- item(a = 0.902, b = c(-1.411, 0.385, 1.79), model = "GRM")
#' plot(ip)
#' plot(ip, category_names = c("Strongly Disagree", "Disagree", "Agree",
#'                             "Strongly Agree"))
#' ip <- item(a = 0.8, b = 1, model = "GRM")
#' plot(ip, category_names = c("Incorrect", "Correct"), legend_title = "Response")
#'
#' # # Change the y-axis label
#' # plot(ip, suppress_plot = TRUE) + ylab("New Label")
#'
plot.Item <- function(x, theta_range = c(-4,4), title = "",
                      suppress_plot = FALSE, category_names = FALSE,
                      legend_title = NULL, ...) {
  theta <- seq(from=theta_range[1], to=theta_range[2], length.out = 501)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if (x@model == 'GRM') {
      icc <- prob(x, theta = theta)
      if (is.null(category_names) ||
          # if the category_names are FALSE.
          (all(is.logical(category_names)) && !all(category_names))) {
        category_labels <- colnames(icc)
      } else category_labels <- category_names
      if (title == "") title <- "Category Response Functions"
      icc <- reshape(data.frame(icc), varying = list(colnames(icc)),
                     v.names = 'p', ids = theta, idvar = 'theta',
                     timevar ='Category', times = category_labels,
                     direction = 'long')
      # Order categories so that the order correctly captured in graph legend.
      icc$Category <- factor(icc$Category, levels = category_labels,
                             ordered = TRUE)
      if (is.null(legend_title)) legend_title <- "Category"
      p <- ggplot(data = icc,
                  aes_string(x = 'theta', y = 'p', color = 'Category')) +
        # Do not show legend if the category_names is NULL
        geom_line(..., show.legend = !is.null(category_names)) +
        xlab(expression("Theta ("*theta*")")) +
        ylab("Probability of Response") +
        ggtitle(title) +
        theme(text = element_text(size=18)) +
        scale_x_continuous(breaks = seq(from = ceiling(theta_range[1]),
                                        to = floor(theta_range[2]), 1)) +
        ylim(0, 1) +
        scale_color_discrete(name=legend_title) +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4))) +
        theme_bw()
    } else {
      ylabel <- "Probability of Correct Response"
      icc <- data.frame(theta = theta, p = prob(ip = x, theta = theta))
      if (title == "") title <- "Item Characteristic Curve"
      # If there is only one item do not print out the legend
      p <- ggplot(data = icc, aes(x = theta, y = p)) +
        geom_line(...) +
        xlab(expression("Theta ("*theta*")")) +
        ylab(ylabel) +
        ylim(0, 1) +
        ggtitle(title) +
        scale_x_continuous(breaks = seq(from = ceiling(theta_range[1]),
                                                 to = floor(theta_range[2]), 1)) +
        theme(text=element_text(size=18)) + theme_bw()
    }
    if (suppress_plot) return(p) else print(p)
  }
}

############################################################################%###
############################# plot.Itempool ###############################%###
############################################################################%###
#' Plot Item Characteristic Curves or Test Characteristic Curve of an
#' \code{Itempool} object
#'
#' @description
#' \code{plot.Itempool} plots the item characteristic curves (item response
#' curves) or test characteristic curve of an \code{\link{Itempool-class}} object.
#'
#' @param x An \code{\link{Itempool-class}} object.
#' @param theta_range The boundaries of x axis.
#' @param tcc If \code{TRUE} a test characteristic curve will be plotted.
#' @param tcc_prop_corr If \code{TRUE}, test characteristic curve will be
#'          show the proportion correct of the test (i.e. the range of y-axis
#'          will be 0-1 instead of 0 to the number of items).
#' @param title Title of the plot. Default is \code{NULL}. If \code{tcc} is
#'          \code{TRUE} it will be 'Test Characteristic Curve', if \code{FALSE}
#'          it will be 'Item Characteristic Curve'.
#' @param suppress_plot If \code{FALSE} the function will print the plot. If
#'          \code{TRUE}, function will return the plot object. Default value is
#'          \code{FALSE}.
#' @param legend_title The title of the plot's legend.
#' @param ... Additional arguments that will be passed to \code{geom_line}
#'
#' @return Depending on the value of \code{suppress_plot} function either prints
#' the item characteristic curve or returns the plot object.
#'
#' @export
#' @importFrom ggplot2 aes aes_string element_text geom_line ggplot ggtitle
#'             guide_legend guides theme theme_bw xlab ylab scale_x_continuous
#'             scale_color_discrete guide_axis
#'
#' @author Emre Gonulates
#'
#' @examples
#' n <- sample(10:15,1)
#' ip <- itempool(a = runif(n, .5, 2), b = rnorm(n), c = runif(n, 0, .3), D = 1)
#' plot(ip)
#' # Additional arguments will passed to geom_line
#' plot(ip, size = .25, alpha = 0.3)
#' # Test Characteristic Curve
#' plot(ip, tcc = TRUE)
#' # Proportion correct for test characteristic curve
#' plot(ip, tcc = TRUE, tcc_prop_corr = TRUE)
#' # # Remove the legend altogether
#' # plot(ip, suppress_plot = TRUE) + theme(legend.position="none")
#' # # Change the labels:
#' # plot(ip, suppress_plot = TRUE) + ylab("Probability") + xlab("Ability Score")
plot.Itempool <- function(x, theta_range = c(-4,4), tcc =FALSE,
                          tcc_prop_corr = FALSE,
                          title = "", suppress_plot = FALSE,
                          legend_title = NULL, ...) {
  theta <- seq(from=theta_range[1], to=theta_range[2], length.out = 501)
  if (title == "")
    title <- ifelse(tcc, "Test Characteristic Curve", "Item Characteristic Curve")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if (tcc) {
      ylabel <- "True Score"
      icc <- data.frame(theta = theta, p = rowSums(prob(ip = x, theta = theta)))
      if (tcc_prop_corr) {
        icc$p <- icc$p/length(x)
        ylabel <- "Proportion Correct"
      }
    } else {
      ylabel <- "Probability of Correct Response"
      icc <- data.frame(itemID = character(0), theta = numeric(0), p = numeric(0))
      for (i in 1:length(x@item_list))
        icc <- rbind(icc, data.frame(itemID = x@item_list[[i]]@id,
                                     theta = theta,
                                     p = prob(ip = x@item_list[[i]], theta = theta)))
    }
    if (is.null(legend_title)) legend_title <- "Item id"
    # If there is only one item do not print out the legend
    if (tcc) {
      p <- ggplot(data = icc, aes_string(x = "theta", y = "p"))
    } else
      p <- ggplot(data = icc, aes_string(x = "theta", y = "p", color = "itemID"))
    p <- p + geom_line(...) +
      xlab(expression("Theta ("*theta*")")) +
      ylab(ylabel) + ggtitle(title) +
      scale_x_continuous(breaks = seq(from = ceiling(theta_range[1]),
                                               to = floor(theta_range[2]), 1)) +
      scale_color_discrete(name=legend_title) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4)))
    p <- p + theme(text=element_text(size=18)) + theme_bw()
    if (suppress_plot) return(p) else print(p)
  }
}

############################################################################%###
############################# plot_info ####################################%###
############################################################################%###
#' Plot Item Information Function
#' @description
#' \code{plot_info} Plots the item information function.
#'
#' @param ip An \code{\link{Item-class}} or \code{\link{Itempool-class}}
#'   object.
#' @param tif If \code{TRUE} a test information plot will be plotted. The
#'          default value is \code{FALSE}.
#' @param theta_range The boundaries of x axis.
#' @param title Title of the plot
#' @param suppress_plot If \code{FALSE} the function will print the plot. If
#'          \code{TRUE}, function will return the plot object. Default value is
#'          \code{FALSE}.
#' @param ... Extra parameters that will pass to \code{geom_line}.
#'
#' @return Depending on the value of \code{suppress_plot} function either prints
#' the item information function or returns the plot object.
#'
#' @export
#' @importFrom ggplot2 aes aes_string element_text geom_line ggplot ggtitle
#'             guide_legend guides theme theme_bw xlab ylab scale_x_continuous
#'             scale_color_discrete
#' @author Emre Gonulates
#'
#' @examples
#' # Plot the information function of an item
#' plot_info(item(b = 1))
#'
#' # Plot information function(s) of an Itempool object
#' n <- sample(10:20,1)
#' ip <- itempool(data.frame(a = runif(n, .5, 2), b = rnorm(n),
#'                              c = runif(n, 0, .3), D = 1))
#' plot_info(ip)
#' plot_info(ip, tif = TRUE)
plot_info <- function(ip, tif = FALSE, theta_range = c(-5,5),
                      title = "", suppress_plot = FALSE,
                      ...)
{
  # Convert ip to Itempool object
  if (!is(ip, "Itempool"))
    tryCatch(
      ip <- itempool(ip),
      error = function(cond) {
        message("\nip cannot be converted to an 'Itempool' object. \n")
        stop(cond)
      })
  if (title == "") title = ifelse(tif, "Test Information Function",
                                  "Item Information Function")
  theta <- seq(from=theta_range[1], to=theta_range[2], length.out = 300)
  info_matrix <- info(ip, theta, tif = tif)
  if (tif) info_data <- data.frame(theta = theta, info = info_matrix) else
    info_data <- data.frame(itemID = rep(ip$id, each = length(theta)),
                            theta = rep(theta, length(ip)),
                            info = as.vector(info_matrix))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if (tif | length(ip) == 1) {
      p <- ggplot(data = info_data, aes_string(x = 'theta', y = 'info'))
    } else
      p <- ggplot(data = info_data, aes_string(x = 'theta', y = 'info', color = 'itemID'))
    p <- p +  geom_line(...) +
      xlab(expression("Theta ("*theta*")")) +
      ylab("Information") + ggtitle(title) +
      theme(text=element_text(size=18))
    if (!tif)
      p <- p + scale_color_discrete(name="Item id") +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4)))
    p <- p + theme_bw()
    if (suppress_plot) return(p) else print(p)
  } else NULL
}

############################################################################%###
############################# plot_resp_loglik #############################%###
############################################################################%###

#' Plot the Log-Likelihood of a response string
#' @description
#' \code{plot_resp_loglik} plots the log-likelihood of a response string.
#'
#' @param ip An \code{\link{Itempool-class}} class object.
#' @param resp The response string
#' @param theta_range The boundaries of x axis.
#' @param title Title of the Plot
#' @param likelihood If \code{TRUE}, likelihood function will be plotted
#'          instead of log-likelihood graph. Default value is \code{FALSE}.
#' @param show_estimate If \code{TRUE} the maximum likelihood ability estimate
#'          will be shown. The default value is \code{TRUE}.
#' @param suppress_plot If \code{FALSE} the function will print the plot. If
#'          \code{TRUE}, function will return the plot object. Default value is
#'          \code{FALSE}.
#' @param text_size The overall text size of the axis and titles. The default
#'          value is 12.
#' @param ... Additional arguments passed to annotate.
#'
#' @return Depending on the value of \code{suppress_plot} function either prints
#' the Log-likelihood function of the response string or returns the plot
#' object.
#'
#' @export
#' @importFrom ggplot2 aes aes_string annotate element_text geom_line
#'             geom_vline ggplot ggtitle theme theme_bw xlab ylab
#' @author Emre Gonulates
#'
#' @section To-do:
#' \itemize{
#' \item Make it to plot multiple test information functions. You can input a
#' list each of which contains item parameters. And the name of the test also.
#' }
#' @examples
#' n <- sample(10:50,1)
#' ip <- itempool(data.frame(a = runif(n, .5, 2), b = rnorm(n),
#'                              c = runif(n, 0, .3), D = 1.7))
#' resp <- sim_resp(ip = ip, theta = rnorm(1))
#' plot_resp_loglik(ip, resp)
#' plot_resp_loglik(ip, resp, text_size = 9)
#' # Format the text of the MLE estimate
#' plot_resp_loglik(ip, resp, size = 3, color = 'blue')
#' # Suppress the MLE estimate
#' plot_resp_loglik(ip, resp, show_estimate = FALSE)
plot_resp_loglik <- function(ip, resp, theta_range = c(-5,5), title = "",
                             likelihood = FALSE, show_estimate = TRUE,
                             suppress_plot = FALSE, text_size = 12, ...) {
  # Convert ip to Itempool object
  if (!is(ip, "Itempool"))
    tryCatch(
      ip <- itempool(ip),
      error = function(cond) {
        message("\nip cannot be converted to an 'Itempool' object. \n")
        stop(cond)
      })
  if (title == "")
    title <- ifelse(likelihood, "Likelihood Graph", "Log-Likelihood Graph")
  theta <- seq(from=theta_range[1], to=theta_range[2], length.out = 201)
  # Prepare a dataset for plot
  value <- sapply(theta, FUN = function(x) resp_loglik(ip = ip, resp = resp,
                                                       theta = x))
  graph_data <- data.frame(theta = theta)
  if (likelihood) graph_data$value <- exp(value) else graph_data$value <- value
  xIntercept <- graph_data$theta[which.max(graph_data$value)]
  yIntercept <- min(graph_data$value)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot(graph_data, aes_string(x = 'theta', y = 'value')) +
      geom_line(size = 1, color = "red") +
      xlab(expression("Theta ("*theta*")")) +
      ylab(ifelse(likelihood, "Likelihood", "Log-Likelihood")) +
      ggtitle(title)
    if (show_estimate)
      p <- p + geom_vline(xintercept = xIntercept, linetype = "dashed") +
        annotate("text", x = xIntercept, y = yIntercept,
                 label = paste("hat(theta) == ", xIntercept),
                 hjust = ifelse(xIntercept > stats::median(theta_range), 1.1, -.1),
                 parse = T, ...)
    p <- p + theme_bw() + theme(text=element_text(size=text_size))
    if (suppress_plot) return(p) else print(p)
  }
}


############################################################################%###
############################# plot_empirical_icc ###########################%###
############################################################################%###
#' Plot Empirical Item or Test characteristic curve
#' @description
#' \code{plot_emprical_icc} plots empirical item or test characteristic curve.
#'
#' @param resp Response matrix.
#' @param item The column number, column name or the 'id' of the  the item that
#'   should be plotted.
#' @param type The type of the graph that will be plotted.
#'   \describe{
#'     \item{\strong{\code{"eicc"}}}{Plot empirical item characteristic curve.
#'       Examinees will be put into bins based on their total raw scores and the
#'       proportion of examinees who correctly answered an item for each bin
#'       will be plotted.}
#'     \item{\strong{"oep"}}{Plot Observed p-values vs. expected p-values
#'       grouped into bins based on total raw scores or theta scores.
#'       This plot requires an \code{\link{Itempool-class}} object.
#'       Optionally, provide \code{theta} vector, otherwise examinee abilities
#'       will be estimated by \code{est_ability(..., type = "eap")}. This will
#'       slow down the plotting function.
#'     }
#'   }
#' @param bins An integer larger than 2 representing of ability groups examinees
#'   should be grouped into. The default is \code{10}. The maximum value of
#'   \code{bins +  1} is the number of possible total scores.
#' @param ip An \code{\link{Itempool-class}} object that is needed for some
#'   plots.
#' @param theta A vector of examinee abilities.
#' @param title Title of the plot
#' @param n_dodge The number of lines the x-axis tick labels should be written
#'   to. This is especially useful if the x-axis tick labels overlap with each
#'   other. The default value is \code{1}, which means all of the labels are
#'   written on the same line.
#' @param suppress_plot If \code{FALSE} the function will print the plot. If
#'          \code{TRUE}, function will return the plot object. Default value is
#'          \code{FALSE}.
#' @param x_axis_scale Set the scale of the x-axis. The default value is
#'   \code{NULL}. For total score it will be defaulted to \code{"percent"}.
#'   \describe{
#'     \item{\strong{\code{"percent"}}}{Percent interval.}
#'     \item{\strong{\code{"number"}}}{Numbers between 1 and \code{bins}}
#'     \item{\strong{\code{"theta"}}}{Theta values equally divided into bins.
#'     the middle value of the bin is shown in the x-axis. For example, if
#'     \code{bins = 10}, the first tick of the x-axis will be the mean of
#'     minimum theta value and tenth percentile theta value.
#'
#'     This is the only option for \code{type = "oep"}. }
#'     }
#'
#' @param ... Extra parameters that will pass to \code{geom_line}.
#'
#' @return Depending on the value of \code{suppress_plot} function either prints
#' the empirical item or test characteristic curve or returns the plot object.
#'
#' @export
#'
#' @importFrom ggplot2 aes aes_string element_text geom_line ggplot ggtitle
#'             guide_legend guides theme theme_bw xlab ylab xlim geom_point
#'             element_blank scale_x_discrete
#' @author Emre Gonulates
#'
#' @examples
#' # Plot the information function of an item
#' resp <- sim_resp(ip = generate_ip(model = "3PL", n = 20),
#'                  theta = rnorm(10000))
#' plot_empirical_icc(resp, 3)
#' # Change the number of bins
#' plot_empirical_icc(resp, 4, bins = 15)
#'
plot_empirical_icc <- function(resp, item, type = "eicc", bins = 10, ip = NULL,
                              theta = NULL, title = "", suppress_plot = FALSE,
                              x_axis_scale = NULL, n_dodge = 1, ...) {
  ########## Check Arguments ########@###
  # resp should be a matrix or data.frame
  if (!inherits(resp, c("matrix", "data.frame")))
    stop("'resp' should be a matrix or a data.frame object.")
  # There should be at least two bins and bins should be a number.
  if (!is.numeric(bins) || as.integer(bins) != bins || bins < 2)
    stop("'bins' should be an integer larger than 2.")
  # Check 'item' argument if it is character it should be an item id, if it
  # is a number it should be a
  item_id <- NULL
  if (is.numeric(item) && (length(item) == 1) && (item <= ncol(resp))) {
    item_col_no <- item
    if (!is.null(colnames(resp))) {
      item_id <- colnames(resp)[item]
    } else if (!is.null(ip) && item %in% ip$id) {
      item_id <- ip$id[item]
    }
  } else if (is.character(item) && (length(item) == 1) &&
            (item %in% colnames(resp) || (!is.null(ip) && item %in% ip$id))) {
    if (item %in% colnames(resp)) {
      item_col_no <- which(item == colnames(resp))
    } else {
      item_col_no <- which(item == colnames(resp))
    }
  } else
    stop("Invalid 'item' argument. Please provide a valid column number, ",
         "column name or the item 'id'. Only provide one 'item'.")
  ### Itempool checks
  # Check whether graph type requires a supported item pool:
  if (type %in% c("oep", "oept")) {
    supported_models <- c(names(Pmodels)[
      sapply(Pmodels, function(x) x$model_family == "UIRT")],
      "GRM", "GPCM", "GPCM2", "PCM")
    if (!is(ip, "Itempool")) {
      stop(paste0("The plot type '", type,
                  "' requires a valid Itempool object."))
      # Currently only dichotomous IRT models supported.
    } else if (!all(ip$model %in% supported_models))  {
      stop(paste0("Currently, only following models are supported by this ",
                  "function: ",
                  paste0("'", supported_models, "'", collapse = ", ")))
    }
  }
  ### x_axis_scale
  if (is.null(x_axis_scale)) x_axis_scale <- "percent"
  if (x_axis_scale == "theta" && !is(ip, "Itempool"))
    stop("When x_axis_scale is 'theta', a valid item pool ('ip') should be ",
         "provided.")
  if (type == "oep" && x_axis_scale != "theta") {
    x_axis_scale <- "theta"
    warning("When type = 'oep', the only available x_axis_scale value is ",
            "'theta'.")
  }
  ### theta checks
  # Check whether graph type requires a theta:
  if (type %in% c("oep", "oept") && !is.numeric(theta)) {
    # Check whether theta can be estimated
    if (is(ip, "Itempool")) {
      theta <- est_ability(ip = ip, resp = resp, method = "eap")$est
    } else
      stop(paste0("The plot type '", type, "' requires a valid theta or an ",
                  "Itempool object to estimate theta."))
  }
  # "eicc" will be calculated using theta estimates instead of total raw scores
  # if 'ip' is available.
  if (type == "eicc" && is(ip, "Itempool")) {
    if (!is.numeric(theta)) {
      theta <- est_ability(ip = ip, resp = resp, method = "eap")$est
    }
  }
  # Theta length should be equal to the number of rows of the resp.
  if (!is.null(theta) && length(theta) != nrow(resp))
    stop("The length of theta should be equal to the number of columns of ",
         "'resp'.")
  # Set the graph title:
  if (title == "")
    title <- paste0("Trace Graph for ", ifelse(
      is.null(colnames(resp)), paste0("Item ", item_col_no),
              colnames(resp)[item_col_no]))
  p <- NULL # by default plot nothing.
  # Convert to data.frame without changing the column names
  temp <- data.frame(resp, check.names = FALSE)
  # Add total score
  if (is.null(theta)) {
    temp$ts <- rowSums(temp, na.rm = TRUE)
  } else temp$ts <- theta
  # Sort the data using the total scores and extract only the relevant column
  # and the total score
  temp <- temp[order(temp$ts), c(item_col_no, ncol(temp))]
  # Remove rows where response is NA
  temp <- temp[!is.na(temp[, 1]), ]

  if (type %in% c("eicc", "oep")) {
    repeat {
      q <- stats::quantile(
        temp$ts, probs = seq(0, 1, length.out = bins+1))[-c(1, bins+1)]
      if (any(duplicated(q)) || any(range(temp$ts, na.rm = TRUE) %in% q)) {
        bins <- bins - 1
        next
      }
      break
    }
    if (x_axis_scale == "percent") {
      names(q) <- paste0(round(as.numeric(gsub("[%]", "", names(q)))), "%")
      temp_labels <- paste0(c("0", names(q)), " - ", c(names(q), "100%"))
    } else if (x_axis_scale == "number") {
      temp_labels <- paste0(1:bins)
    } else if (x_axis_scale == "theta") {
      temp_labels <- unname(round((c(q, max(temp$ts)) +
                                     c(min(temp$ts), q)) / 2, 2))
    } else stop("Invalid 'x_axis_scale' value.")

    q <- c(min(temp$ts, na.rm = TRUE), q, max(temp$ts, na.rm = TRUE))
    temp$bin <- cut(temp$ts, breaks = q, include.lowest = TRUE,
                    labels = temp_labels)
    temp[ , 1] <- factor(temp[, 1], ordered = TRUE)

    temp <- as.data.frame(t(prop.table(table(Response = temp[, 1],
                                             bin = temp$bin), margin = 2)))

    # If there are just two categories, remove the first category.
    if (length(unique(temp$Response)) == 2)
        temp <- temp[temp$Response == unique(temp$Response)[2], ]

    # Set x-axis label
    x_label <- ifelse(is.null(theta), "Ability Group (Raw Score)",
                      "Ability Group (Theta Score)")
    if (type == "eicc") {
      if (length(unique(temp$Response)) > 1) {
      # if (length(empty_element) > 2) {
      #   temp <- stats::reshape(
      #     data = temp, varying = list(names(empty_element)),
      #     # data = temp, varying = list(names(temp)[-ncol(temp)]),
      #     v.names = "Proportion", direction = "long", timevar = "Response",
      #     times = names(temp1[[1]]))
      #   p <- ggplot(temp, aes_string(x = "bin", y = "Proportion",
      #                                color = "Response")) +
      #     geom_point(size = 2) +
      #     geom_line(aes_string(group = "Response"), ...)
        p <- ggplot(temp, aes_string(x = "bin", y = "Freq",
                                     color = "Response")) +
          geom_point(size = 2) +
          geom_line(aes_string(group = "Response"), ...)
      } else {
        # colnames(temp) <- c("Proportion", "n", "bin")
        p <- ggplot(temp, aes_string(x = "bin", y = "Freq", group = 1)) +
          geom_point(size = 2) + geom_line(...)
      }
    } else if (type == "oep") {
      colnames(temp)[colnames(temp) == "Freq"] <- "Observed"
      temp$Expected <- resp_lik(
        ip = ip[[item_col_no]],
        resp = as.integer(as.character(temp$Response)),
        theta = as.numeric(as.character(temp$bin)))
      temp <- stats::reshape(
        data = temp, varying = c("Observed", "Expected"), v.names = "pvalue",
        direction = "long", timevar = "Type", times = c("Observed", "Expected"))
      if (x_axis_scale == "theta")
        temp$bin <- as.numeric(as.character(temp$bin))
      if (length(unique(temp$Response)) > 2) {
        p <- ggplot(temp, aes_string(
          x = "bin", y = "pvalue", color = "Response", linetype = "Type",
          group = "interaction(Type, Response)")) +
          geom_point() + geom_line(...)
      } else {
        p <- ggplot(temp, aes_string(x = "bin", y = "pvalue", color = "Type")) +
          geom_point(group = "Type", size = 2) +
          geom_line(aes_string(group = "Type"), ...)
      }
    }
  } else stop("This graph type has not been implemented yet.")
  # Add common graph elements
  if (!is.null(p)) {
    p <- p + xlab(x_label) + ylab("Proportion Correct") + ggtitle(title) +
      scale_x_discrete(guide = guide_axis(n.dodge = n_dodge)) +
      ylim(c(0, 1)) + theme_bw()
    if (type == "oep") p <- p + theme(legend.title=element_blank())
  }
  if (suppress_plot) return(p) else print(p)
}


############################################################################%###
############################# plot_distractor_icc ##########################%###
############################################################################%###
#' Plot Empirical Item or Test characteristic curve
#' @description
#' \code{plot_empirical_icc} plots empirical item or test characteristic curve.
#'
#' @param raw_resp Raw response matrix.
#' @param item The column number, column name or the 'id' of the  the item that
#'   should be plotted.
#' @param key A vector of answer key.
#' @param bins An integer larger than 2 representing of ability groups examinees
#'   should be grouped into. The default is \code{10}. The maximum value of
#'   \code{bins +  1} is the number of possible total scores.
#' @param ip An \code{\link{Itempool-class}} object that is needed for some
#'   plots. If \code{ip} provided and \code{theta} is not provided, then
#'   ability will be estimated using EAP method with prior mean 0 and prior
#'   standard deviation of 1. This is a slower method depending on the size of
#'   the data.
#' @param theta A vector of examinee abilities. If \code{theta} values provided
#'   the bins are formed using them instead of sum scores.
#' @param x_axis_scale Set the scale of the x-axis. The default value is
#'   \code{NULL}. For if sum score is used scale will be defaulted to
#'   \code{"percent"}, Otherwise if valid \code{theta} or \code{ip} arguments
#'     provided  the scale defaults to \code{"theta"}.
#'   \describe{
#'     \item{\strong{\code{"percent"}}}{Percent interval.}
#'     \item{\strong{\code{"number"}}}{Numbers between 1 and \code{bins}.}
#'     \item{\strong{\code{"theta"}}}{Theta values equally divided into bins.
#'     the middle value of the bin is shown in the x-axis. For example, if
#'     \code{bins = 10}, the first tick of the x-axis will be the mean of
#'     minimum theta value and tenth percentile theta value.}
#'     }
#' @param title Title of the plot
#' @param n_dodge The number of lines the x-axis tick labels should be written
#'   to. This is especially useful if the x-axis tick labels overlap with each
#'   other. The default value is \code{1}, which means all of the labels are
#'   written on the same line.
#' @param suppress_plot If \code{FALSE} the function will print the plot. If
#'          \code{TRUE}, function will return the plot object. Default value is
#'          \code{FALSE}.
#' @param ... Extra parameters that will pass to \code{geom_line}.
#'
#' @return Depending on the value of \code{suppress_plot} function either prints
#' the proportion of examinees in each bin respond to each distractor  or
#' returns the plot object.
#'
#' @export
#' @importFrom ggplot2 aes aes_string element_text geom_line ggplot ggtitle
#'             guide_legend guides theme theme_bw xlab ylab geom_point
#'             scale_linetype_discrete scale_x_discrete guide_axis
#' @author Emre Gonulates
#'
#' @examples
#' n_item <- 10 # sample(8:12, 1)
#' n_theta <- 10000 # sample(100:200, 1)
#' raw_resp <- matrix(sample(LETTERS[1:4], n_item * n_theta, replace = TRUE),
#'                    nrow = n_theta, ncol = n_item,
#'                    dimnames = list(paste0("Examinee-", 1:n_theta),
#'                                    paste0("Item-", 1:n_item)))
#' key <- sample(LETTERS[1:4], n_item, replace = TRUE)
#' plot_distractor_icc(raw_resp, 3, key)
#' # Change the number of bins
#' plot_distractor_icc(raw_resp, 3, key, bins = 15)
#'
plot_distractor_icc <- function(raw_resp, item, key, bins = 10, ip = NULL,
                                theta = NULL, x_axis_scale = NULL,
                                title = "", n_dodge = 1, suppress_plot = FALSE,
                                ...) {
  ########## Check Arguments ########@###
  # raw_resp should be a matrix or data.frame
  if (!inherits(raw_resp, c("matrix", "data.frame")))
    stop("'raw_resp' should be a matrix or a data.frame object.")
  # There should be at least two bins and bins should be a number.
  if (!is.numeric(bins) || as.integer(bins) != bins || bins < 2)
    stop("'bins' should be an integer larger than 2.")
  # Check 'item' argument if it is character it should be an item id, if it
  # is a number it should be a
  item_id <- NULL
  if (is.numeric(item) && (length(item) == 1) && (item <= ncol(raw_resp))) {
    item_col_no <- item
    if (!is.null(colnames(raw_resp))) item_id <- colnames(raw_resp)[item]
  } else if (is.character(item) && (length(item) == 1) &&
             item %in% colnames(raw_resp)) {
    if (item %in% colnames(raw_resp)) {
      item_col_no <- which(item == colnames(raw_resp))
    } else {
      item_col_no <- which(item == colnames(raw_resp))
    }
  } else
    stop("Invalid 'item' argument. Please provide a valid column number, ",
         "column name or the item 'id'. Only provide one 'item'.")
  # Valid key
  if (length(key) != ncol(raw_resp))
    stop("Invalid 'key'. The length of 'key' should be the same as the number ",
         "of columns of 'raw_resp'.")

  # This conversion is necessary to implement a subsetting method that is
  # valid for both matrix, data.frame and tibbles.
  if (inherits(raw_resp, "matrix"))
    raw_resp <- data.frame(raw_resp, check.names = FALSE)
  # Set the graph title:
  if (title == "")
    title <- paste0("Trace Lines for ", ifelse(
      is.null(colnames(raw_resp)), paste0("Item ", item_col_no),
              colnames(raw_resp)[item_col_no]))

  # Sort the data using the total scores and extract only the relevant column
  # and the total score
  # Check theta and Itempool
  if (is.numeric(theta) && length(theta) == nrow(raw_resp)) {
    temp <- data.frame(ts = theta)
  } else if (is(ip, "Itempool")) {
    temp <- data.frame(ts = est_ability(
      ip = ip, resp = score_raw_resp(raw_resp, key), method = "eap"))
  } else { # Calculate the sum scores
    temp <- data.frame(ts = rowSums(score_raw_resp(raw_resp, key),
                                    na.rm = TRUE))
  }
  temp$response <- ifelse(raw_resp[[item_col_no]] == key[item_col_no],
                          paste0(key[item_col_no], "*"),
                          raw_resp[[item_col_no]])
  # Remove rows where response is NA
  temp <- temp[!is.na(temp$response), ]
  # Set the default x_axis_scale:
  if (is.null(x_axis_scale)) {
    if (is.null(theta) && is.null(ip)) x_axis_scale <- "percent" else
      x_axis_scale <- "theta"
  } else if (x_axis_scale == "theta" && (is.null(theta) && is.null(ip)) ) {
    message("x_axis_scale = 'theta' is available only with valid 'theta' or ",
            "'ip'.")
    x_axis_scale <- "percent"
  } else if (!x_axis_scale %in% c("percent", "theta", "number")) {
    stop("Invalid 'x_axis_scale' value.")
  }
  # Set the x-axis label
  if (x_axis_scale == "theta" || !is.null(theta) || !is.null(ip)) {
    x_label <- "Ability Group (Theta Score)"
  } else x_label <- "Ability Group (Sum Score)"

  # Decrease bin number until there are no empty bins
  repeat {
    q <- stats::quantile(temp$ts,
                         probs = seq(0, 1, length.out = bins+1))[-c(1, bins+1)]
    if (any(duplicated(q)) || any(range(temp$ts, na.rm = TRUE) %in% q)) {
      bins <- bins - 1
      next
    }
    break
    # if (length(unique(temp$bin)) == bins || bins == 2) break else
    #   bins <- bins - 1
  }

  if (x_axis_scale == "percent") {
    names(q) <- paste0(round(as.numeric(gsub("[%]", "", names(q)))), "%")
    temp_labels <- paste0(c("0", names(q)), " - ", c(names(q), "100%"))
  } else if (x_axis_scale == "number") {
    temp_labels <- paste0(1:bins)
  } else if (x_axis_scale == "theta") {
    temp_labels <- unname(round((c(q, max(temp$ts)) +
                                   c(min(temp$ts), q)) / 2, 2))
  } else stop("Invalid 'x_axis_scale' value.")

  q <- c(min(temp$ts, na.rm = TRUE), q, max(temp$ts, na.rm = TRUE))

  temp$bin <- cut(temp$ts, breaks = q, include.lowest = TRUE,
                  labels = temp_labels)
  temp <- as.data.frame(t(prop.table(table(Response = temp$response,
                                           bin = temp$bin), margin = 2)))
  temp$incorrect <- !grepl("[*]$", temp$Response)

  p <- ggplot(temp, aes_string(x = "bin", y = "Freq",
                                color = "Response")) +
    geom_point(alpha = 0.75, size = 2) +
    # geom_point(aes_string(size = "n"), alpha = 0.75) +
    geom_line(aes_string(group = "Response", linetype = "incorrect"),
              size = 1, alpha = 0.75, ...) +
    scale_linetype_discrete(guide = FALSE) +
    xlab(x_label) + ylab("Proportion") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    ggtitle(title) + ylim(c(0, 1)) + theme_bw()
  if (suppress_plot) return(p) else print(p)
}

#' ############################################################################%###
#' ############################# hist.Itempool ###############################%###
#' ############################################################################%###
#' #' Plot Test Information Function
#' #' @description
#' #' \code{hist.Itempool} Plots the histogram of the item pool.
#' #'
#' #' @param ip An \code{\link{Itempool-class}} object.
#' #' @param parameter Which item pool parameter to plot.
#' #'                  The choices are: 'a': item discrimination parameter, $a$.
#' #'                                   'b': item difficulty parameter, $b$.
#' #'                                   'c': pseudo-guessing parameter, $c$.
#' #'                  The default value is 'b'.
#' #' @param D Scaling constant.
#' #' @param xlim If
#' #' @param separateContent Whether to separate the content areas in the
#' #'                        histogram. If there are no content areas in the
#' #'                        item pool, there won't be any separation.
#' #'                        Default is \code{FALSE}.
#' #' @param alpha The transparency of the bins. It takes values between 0 and 1.
#' #'              By default, if \code{separateContent} is \code{TRUE} then
#' #'              the value will be 0.5, it will be 1 otherwise.
#' #' @param addTIFline If \code{TRUE}, a test information line will be added
#' #'                   to the histogram. There will be secondary y-axis to the
#' #'                   right of the graph. If \code{separateContent} is
#' #'                   \code{TRUE}, the the TIF will be added to the each
#' #'                   content area. The default is \code{FALSE}.
#' #' @param suppress_plot If \code{FALSE} the function will print the plot. If
#' #'          \code{TRUE}, function will return the plot object. Default value is
#' #'          \code{FALSE}.
#' #' @param ... The arguments that will be passed to the \code{plotHistogram}
#' #'            function.
#' #'
#' #' @return Depending on the value of \code{suppress_plot} function either prints
#' #' a histogram of item pool or returns the plot object.
#' #'
#' #' @export
#' #' @author Emre Gonulates
#' #' @seealso \code{\link{plotHistogram}}
#' #'
#' #' @section To-do:
#' #' \itemize{
#' #'   \item Add title example.
#' #'   \item Add suppress_plot=FALSE example.
#' #' }
#' #' @examples
#' #' \donttest{
#' #' n <- sample(10:50,1);   D = 1
#' #' ip <- itempool(data.frame(a = runif(n, .5, 2), b = rnorm(n),
#' #'                              c = runif(n, 0, .3)))
#' #' plotItemPool(ip)
#' #'
#' #' n <- c(15, 40, 10)
#' #' ip <- itempool(data.frame(
#' #'   a = runif(sum(n), 0.75, 2),
#' #'   b = c(sort(rnorm(n = n[1], mean = 1, sd = 0.3)),
#' #'         sort(rnorm(n = n[2], mean = 0, sd = 0.5)),
#' #'         sort(rnorm(n = n[3], mean = -1, sd = 0.4))),
#' #'   c = runif(sum(n), 0.01, 0.3)
#' #'   ),
#' #'   id = paste0("Item-",1:sum(n)),
#' #'   content = c(rep("Trigonometry", n[1]), rep("Algebra", n[2]),
#' #'               rep("Arithmetic", n[3]))
#' #' )
#' #'
#' #' plotItemPool(ip = ip)
#' #' plotItemPool(ip = ip, binWidth = 0.05)
#' #' plotItemPool(ip = ip, separateContent = TRUE)
#' #' plotItemPool(ip = ip, separateContent = TRUE, addTIFline = TRUE)
#' #' plotItemPool(ip = ip, addTIFline = TRUE)
#' #' plotItemPool(ip = ip, variableLabel = "Item Difficulty Parameter")
#' #' plotItemPool(ip = ip, parameter = 'c', binWidth = 0.025)
#' #' plotItemPool(ip = ip, parameter = 'c', separateContent = TRUE, alpha = .2)
#' #' plotItemPool(ip = ip, parameter = 'a', binWidth = 0.025)
#' #' plotItemPool(ip = ip, parameter = 'a', binWidth = 0.05, separateContent = T)
#' #' plotItemPool(ip = ip, separateContent = T, addTextBox = c("n", 'meanX', 'sdX'))
#' #' plotItemPool(ip = ip, separateContent = TRUE, theme = "theme_bw")
#' #'
#' #' p <- plotItemPool(ip = ip, suppress_plot = TRUE)
#' #'
#' #' # Add another item pool
#' #' ip2 <- itempool(data.frame(a = runif(n[2], .5, 2), b = rnorm(n[2]),
#' #'                               c = runif(n[2], 0, .3)))
#' #' plotItemPool(ip = ip2, binWidth = 0.05)
#' #' plotItemPool(ip = ip2, addTIFline = TRUE, variableLabel = "Item Difficulty")
#' #' plotItemPool(ip = ip2, binWidth = 0.1, seperateContent = TRUE)
#' #' }
#' hist.Itempool <- function(ip, parameter = 'b', xlim = NULL,
#'                           separateContent = FALSE, alpha = NULL,
#'                           addTIFline = FALSE, suppress_plot = FALSE, ...) {
#'   # Convert ip to Itempool object
#'   if (!is(ip, "Itempool"))
#'     tryCatch(
#'       ip <- itempool(ip),
#'       error = function(cond) {
#'         message("\nip cannot be converted to an 'Itempool' object. \n")
#'         stop(cond)
#'       })
#'   argList <- list()
#'   # Determine x-axis label
#'   if (!hasArg(variableLabel)) {
#'     if (parameter == 'b') {variableLabel <- "Item Difficulty (b)"
#'     } else if (parameter == 'a') {variableLabel <- "Item Discrimination (a)"
#'     } else if (parameter == 'c') variableLabel <- "Guessing Parameter (c)"
#'     argList <- c(argList, variableLabel = variableLabel)
#'   }
#'   # Determine bin width
#'   if (!hasArg(binWidth)) {
#'     if (parameter == 'b') {binWidth <- 0.1
#'     } else if (parameter == 'a') {binWidth <- 0.05
#'     } else if (parameter == 'c') binWidth <- 0.02
#'     argList <- c(argList, binWidth = binWidth)
#'   }
#'   ipData <- data.frame(getParameters(object = ip, parameterNames = parameter))
#'   colnames(ipData) <- parameter
#'   # Decide whether to separate histogram by content
#'   separateContent <- separateContent & !is.null(ip$content)
#'   # Add content column if separateContent is TRUE
#'   if (separateContent)
#'     ipData$content <- ip$content
#'   alpha <- ifelse(is.null(alpha), ifelse(separateContent, 0.5, 1), alpha)
#'   # Here, since I'm defining variableLabel and binWidth within this function,
#'   # to prevent any clash, I redefined '...'.
#'   # https://stackoverflow.com/a/17390262/2275286
#'   inargs <- list(...)
#'   argList[names(inargs)] <- inargs
#'   if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     p <- do.call(plotHistogram,
#'                  c(list(data = ipData, variable = parameter, alpha = alpha,
#'                         group = switch(separateContent + 1, NULL, "content"),
#'                         legendTitle = switch(separateContent + 1, NULL, "Content"),
#'                         suppress_plot = FALSE), argList)
#'                  )
#'     if  (parameter == 'b' & addTIFline)
#'     {
#'       if (is.null(xlim))
#'         xlim <- seq(from = -4, to = 4, length.out = 201)
#'       else xlim <- seq(from = xlim[1], to = xlim[2], length.out = 201)
#'       if (separateContent)
#'       {
#'         infoData <- data.frame(theta = numeric(0), information = numeric(0),
#'                                content = character(0))
#'         for (i in sort(unique(ip$content)))
#'           infoData <- rbind(infoData, data.frame(
#'             theta = xlim,
#'             information = info(ip = getContentFast(ip = ip, content = i),
#'                                theta = xlim, tif = TRUE),
#'             content = i))
#'       } else {
#'         infoData <- data.frame(theta = xlim, information = info(
#'           ip = ip, theta = xlim, tif = TRUE))
#'       }
#'       # Find the highest point in histogram
#'       max_hist <- max(ggplot_build(p)$data[[1]]$count)
#'       # Find the highest point in info graph
#'       max_info <- max(infoData$information)
#'       # Find the ratio of max_info/max_hist:
#'       expansion_ratio <- round(max_info/max_hist, 1)
#'       infoData$information <- infoData$information/expansion_ratio
#'
#'       # https://stackoverflow.com/q/3099219/2275286
#'       p <- p + geom_line(data = infoData, aes_string(
#'         x = 'theta', y = 'information',
#'         color = switch(separateContent + 1, NULL, "content")), size = 1) +
#'         scale_y_continuous(name = "Count", sec.axis = sec_axis(
#'           ~.*expansion_ratio, name = "Item Pool Information")) +
#'         scale_color_discrete(name = "Content")
#'     }
#'     if (suppress_plot) return(p) else print(p)
#'   }
#' }
