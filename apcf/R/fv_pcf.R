#' Class fv_pcf: Function Value Table for PCFs
#'
#' Advanced Use Only. This low-level function creates an object of class
#' "fv_pcf" from raw numerical data.
#'
#' @param df A data frame with at least 2 columns named 'r' and 'g' containing
#'        the values of the function argument (r) and the corresponding values
#'        (g). Usually the upper 'upr' and lower 'lwr' bounds of a pointwise
#'        critical envelope are passed along as well.
#' @param n_sim Integer. Number of generated simulated patterns used for
#'        computing the envelope
#' @param n_rank Integer. Rank of the envelope value amongst the n_sim
#'        simulated values. A rank of 1 means that the minimum and maximum
#'        simulated values will be used.
#' @param correc String. Choice of edge correction (eg. "Ripley").
#' @param kernel String. Choice of smoothing kernel (eg. "epanechnikov").
#' @param stoyan Bandwidth coefficient used in smoothing kernel.
#' @param bw Bandwidth used in smoothing kernel.
#' @param x,obj,object an R object, preferably of class `fv_pcf`
#' @param ... additional parameter
#'
#' @return An object of class `fv_pcf`.
#'
#' @export
fv_pcf <- function(df, n_sim, n_rank, correc, kernel, stoyan, bw){
  stopifnot(is.data.frame(df))
  if(n_sim %% 1 != 0)
    stop("n_sim must be an integer")

  if(n_sim < 1){
    stop("n_rank must be >= 1")
  } else if(n_sim == 1){
    if(length(df) < 2 || !all(c("r", "g") %in% names(df)))
      stop("'df' must have at least columns named r, g")
  } else {
    if(length(df) < 4 || !all(c("r", "g", "lwr", "upr") %in% names(df)))
      stop("'df' must have at least columns named r, g, lwr, upr")
  }

  if(n_rank %% 1 != 0)
    stop("n_rank must be an integer")
  if(n_rank < 1){
    stop("n_rank must be >= 1")
  } else if(n_sim  > 1 && !(n_rank < n_sim/2)){
    stop("n_rank mus be < n_sim/2")
  }
  alpha <- (2 * n_rank)/(n_sim+1)

  stopifnot(is.character(correc))
  stopifnot(is.character(kernel))
  stopifnot(is.numeric(stoyan))
  stopifnot(is.numeric(bw))

  attr(df, "n_sim")  <- n_sim
  attr(df, "n_rank") <- n_rank
  attr(df, "alpha")  <- alpha
  attr(df, "correc") <- correc
  attr(df, "kernel") <- kernel
  attr(df, "stoyan") <- stoyan
  attr(df, "bw")     <- bw
  class(df) <- c("fv_pcf", class(df))

  return(df)
}


#' @rdname fv_pcf
#' @export
is.fv_pcf <- function(obj){
  inherits(obj, "fv_pcf")
}


# crudely adopted from spatstat (v1.48-0)
#' @rdname fv_pcf
#' @export
print.fv_pcf <- function(x, ...){
  cat("PCF with pointwise critical envelopes\n")
  cat(paste0("Obtained from ", attr(x, "n_sim"), " simulations"))
  if(attr(x, "n_rank") > 1)
    cat(paste0(" with ", attr(x, "n_rank")-1, " disregarded at top and bottom\n"))
  else
    cat("\n")
  cat(paste0("Edge correction: ", dQuote(attr(x, "correc")), "\n"))
  cat(paste0("Alternative: ", dQuote("two.sided"), "\n"))
  cat(paste0("Significance level of pointwise Monte Carlo test: ",
             attr(x, "alpha"), "\n"))
  cat("\n")
  print.data.frame(x)
}


# crudely adopted from spatstat (v1.48-0)
#' @rdname fv_pcf
#' @export
summary.fv_pcf <- function(object, ...){
  n_rank <- attr(object, "n_rank")

  cat("PCF with pointwise critical envelopes\n")

  cat(paste0("Obtained from ", attr(object, "n_sim"), " simulations"))
  if(n_rank > 1)
    cat(paste0(" with ", n_rank - 1, " disregarded at top and bottom\n"))
  else
    cat("\n")

  cat(paste0("Edge correction: ", dQuote(attr(object, "correc")), "\n"))
  cat(paste0("Alternative: ", dQuote("two.sided"), "\n"))

  lo.ord <- if (n_rank == 1)
    "minimum"
  else
    paste(ordinal(n_rank), "smallest")

  hi.ord <- if (n_rank == 1)
    "maximum"
  else
    paste(ordinal(n_rank), "largest")
  cat("Upper envelope: pointwise", hi.ord, "of simulated curves\n")
  cat("Lower envelope: pointwise", lo.ord, "of simulated curves\n")
  cat(paste0("Significance level of pointwise Monte Carlo test: ",
             attr(object, "alpha"), "\n"))
}


# taken from spatstat (v1.48-0)
ordinal <- function(k){
  last <- abs(k)%%10
  lasttwo <- abs(k)%%100
  isteen <- (lasttwo > 10 & lasttwo < 20)
  ending <- ifelse(isteen, "th",
                   ifelse(last == 1, "st",
                          ifelse(last == 2, "nd",
                                 ifelse(last == 3, "rd",
                                        "th"))))
  return(paste0(k, ending))
}
