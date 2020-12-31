#' Print object of class \code{smerc_cluster}.
#'
#' Print \code{smerc_cluster} object
#'
#' @param x An object of class \code{smerc_cluster}.
#' @inheritDotParams base::print
#' @param extra A logical value. Default is \code{FALSE}.
#' \code{TRUE} indicates that extra information should be
#' printed.
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(x, y))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 49,
#'                 longlat = TRUE, alpha = 0.12)
#' out
print.smerc_cluster = function(x, ..., extra = FALSE) {
  if (length(extra) != 1 | !is.logical(extra)) {
    stop("extra must be a single logical value")
  }
  if (requireNamespace("crayon", quietly = TRUE)) {
    print_smerc_cluster_crayon(x, extra)
  } else {
    cat(paste(("method:"), (x$method)), sep = "\n")
    rel_param = x$rel_param
    if (!is.null(rel_param$type)) {
      cat(paste(("statistic:"),
                    (rel_param$type)), sep = "\n")
    }
    if (!is.null(rel_param$simdist)) {
      cat(paste(("simulation:"),
                    (rel_param$simdist)), sep = "\n")
    }
    if (!is.null(rel_param$nsim)) {
      cat(paste(("realizations:"),
                    (rel_param$nsim)), sep = "\n")
    }
    if (!is.null(rel_param$ubpop)) {
      cat(paste("population upperbound: ",
                rel_param$ubpop * 100, "%", sep = ""),
          sep = "\n")
    }
    if (!is.null(rel_param$ubd)) {
      cat(paste("distance upperbound: ",
                rel_param$ubd * 100, "%", sep = ""),
          sep = "\n")
    }
    if (!is.null(rel_param$min.cases)) {
      cat(paste(("minimum cases:"),
                    (rel_param$min.cases)), sep = "\n")
    }
    if (!is.null(rel_param$a_penalty)) {
      cat(paste(("shape penalty (a):"),
                    (rel_param$a_penalty)), sep = "\n")
    }
    if (!is.null(rel_param$shapes)) {
      cat(paste("shapes:", paste(rel_param$shapes, collapse = ', ')), sep = "\n")
    }
    if (!is.null(rel_param$nangles)) {
      cat(paste("angles per shape:", paste(rel_param$nangles, collapse = ', ')), sep = "\n")
    }
    if (!is.null(rel_param$cstar)) {
      cat(paste(("case radius:"),
                    (rel_param$cstar)), sep = "\n")
    }
    if (!is.null(rel_param$nstar)) {
      cat(paste(("population radius:"),
                    (rel_param$nstar)), sep = "\n")
    }
    if (!is.null(rel_param$modified)) {
      cat(paste(("modified p-value:"),
                    (rel_param$modified)), sep = "\n")
    }
    if (!is.null(rel_param$k)) {
      cat(paste(("number of neighbors:"),
                    (rel_param$k)), sep = "\n")
    }
    if (!is.null(rel_param$alpha1)) {
      cat(paste(("middle p-value threshold:"),
                    (rel_param$alpha1)), sep = "\n")
    }
    if (extra) {
      cat(paste(("significance level:"),
                    (x$alpha)), sep = "\n")
      dtype = ifelse(x$longlat, "great circle", "euclidean")
      cat(paste(("distance:"),
                    (dtype)), sep = "\n")
      cat(paste(("total regions:"),
                    (x$number_of_regions)), sep = "\n")
      cat(paste(("total cases:"),
                    round(x$total_cases, 2)), sep = "\n")
      cat(paste(("total population:"),
                    (x$total_population)), sep = "\n")
      cat(paste(("cases per 100,000 persons:"),
                    round(x$cases_per_100k, 2)), sep = "\n")
    }
  }
}

#' Print smerc_cluster in color
#'
#' Print smerc_cluster in color using crayon package
#'
#' @param x A smerc_cluster object
#' @param extra Should extra information be printed?
#' @return NULL
#' @noRd
print_smerc_cluster_crayon = function(x, extra) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    message(paste(crayon::blue("method:"),
                  crayon::magenta(x$method)))
    rel_param = x$rel_param
    if (!is.null(rel_param$type)) {
      message(paste(crayon::blue("statistic:"),
                    crayon::magenta(rel_param$type)))
    }
    if (!is.null(rel_param$simdist)) {
      message(paste(crayon::blue("simulation:"),
                    crayon::magenta(rel_param$simdist)))
    }
    if (!is.null(rel_param$nsim)) {
      message(paste(crayon::blue("realizations:"),
                    crayon::magenta(rel_param$nsim)))
    }
    if (!is.null(rel_param$ubpop)) {
      message(paste(crayon::blue("population upperbound:"),
                    crayon::magenta(rel_param$ubpop * 100, "%",
                                    sep = "")))
    }
    if (!is.null(rel_param$ubd)) {
      message(paste(crayon::blue("distance upperbound:"),
                    crayon::magenta(rel_param$ubd * 100, "%",
                                    sep = "")))
    }
    if (!is.null(rel_param$min.cases)) {
      message(paste(crayon::blue("minimum cases:"),
                    crayon::magenta(rel_param$min.cases)))
    }
    if (!is.null(rel_param$a_penalty)) {
      message(paste(crayon::blue("shape penalty (a):"),
                    crayon::magenta(rel_param$a_penalty)))
    }
    if (!is.null(rel_param$shapes)) {
      message(paste(crayon::blue("shapes:"),
                    crayon::magenta(paste(rel_param$shapes, collapse = ', '))))
    }
    if (!is.null(rel_param$nangles)) {
      message(paste(crayon::blue("angles per shape:"),
                    crayon::magenta(paste(rel_param$nangles, collapse = ', '))))
    }
    if (!is.null(rel_param$cstar)) {
      message(paste(crayon::blue("case radius:"),
                    crayon::magenta(rel_param$cstar)))
    }
    if (!is.null(rel_param$nstar)) {
      message(paste(crayon::blue("population radius:"),
                    crayon::magenta(rel_param$nstar)))
    }
    if (!is.null(rel_param$modified)) {
      message(paste(crayon::blue("modified p-value:"),
                    crayon::magenta(rel_param$modified)))
    }
    if (!is.null(rel_param$k)) {
      message(paste(crayon::blue("number of neighbors:"),
                    crayon::magenta(rel_param$k)))
    }
    if (!is.null(rel_param$alpha1)) {
      message(paste(crayon::blue("middle p-value threshold:"),
                    crayon::magenta(rel_param$alpha1)))
    }
    if (extra) {
      message(paste(crayon::blue("significance level:"),
                    crayon::magenta(x$alpha)))
      dtype = ifelse(x$longlat, "great circle", "euclidean")
      message(paste(crayon::blue("distance:"),
                    crayon::magenta(dtype)))
      message(paste(crayon::blue("total regions:"),
                    crayon::magenta(x$number_of_regions)))
      message(paste(crayon::blue("total cases:"),
                    crayon::magenta(round(x$total_cases, 2))))
      message(paste(crayon::blue("total population:"),
                    crayon::magenta(x$total_population)))
      message(paste(crayon::blue("cases per 100,000 persons:"),
                    crayon::magenta(round(x$cases_per_100k, 2))))
    }
  }
}


