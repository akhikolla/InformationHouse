
# ==== details ====
details <- function(cond, x, 
  what = c("inus", "cyclic", "exhaustiveness", "faithfulness", "coherence", "redundant"),
  cycle.type = c("factor", "value")){
  cycle.type <- match.arg(cycle.type)
  .det(x, noblanks(cond), what = what, cycle.type = cycle.type)
}


# ==== .det() ====
# Switch order of first 2 args to provide dispatching on x
# Generic function
.det <- function(x, cond, ...) UseMethod(".det")

# ==== Default method for class 'cti' ====
#  Converts x to a configTable and then calls the .det.cti
.det.default <- function(x, cond, 
    what = c("inus", "cyclic", "exhaustiveness", "faithfulness", "coherence", "redundant"),
    type, cycle.type){
  if (!inherits(x, "configTable")){
    if (missing(type)) type <- "cs"
    x <- configTable(x, type = type)
  }
  .det.cti(ctInfo(x), cond = cond, what = what, available = x$details, cycle.type)
}

# ==== Method for class 'cti' ====
#   x        cti
#   cond     character vector with the csf
#   what     char-vector: list of required measures
#   in.csf   special shortcut for the case where .det.cti is called from within csf()
# value: A data.frame with columns according to what
.det.cti <- function(x, cond, 
    what = c("inus", "cyclic", "exhaustiveness", "faithfulness", "coherence", "redundant"),
    available = what, cycle.type, in.csf = FALSE){
  what <- clarify_details(what, available = available)
  qc_ct <- qcond_csf(cond, x$scores, flat = TRUE)
  if (any(c("redundant", "inus", "exhaustiveness", "faithfulness") %in% what)){
    x_full <- full.ct(x)
    qc_full <- qcond_csf(cond, x_full$scores, flat = TRUE)
  }
  out <- data.frame(matrix(numeric(0), length(cond)), row.names = cond)
  if ("redundant" %in% what || (in.csf && "inus" %in% what)){
    red <- m_any(.redund.cti(x_full, cond, simplify = FALSE, full = TRUE, 
                             qc_full = qc_full, names = FALSE))
  }
  if ("inus" %in% what){
    if (in.csf){
      out$inus <- !red
      const <- rep(FALSE, length(cond))
      const[out$inus] <- apply(qcond_csf(cond[out$inus], x_full$scores, force.bool = TRUE), 2, isConstant)
      out$inus[const] <- FALSE
      constFact <- constFact(cond[out$inus], x_full)
      out$inus[out$inus][constFact] <- FALSE
      partStrRed <- partiallyRedundant(extract_asf(cond[out$inus]), x_full)
      out$inus[out$inus][partStrRed] <- FALSE
    } else {
      out$inus <- .inus.cti(x_full, cond, full = TRUE)
    }
  }
  if ("cyclic" %in% what){
    out$cyclic <- cyclic(cond, cycle.type = cycle.type, use.names = FALSE)
  }
  ex_ff <- c("exhaustiveness", "faithfulness") %in% what
  if (any(ex_ff)){
    exf <- exff.cti(x, cond, names = FALSE, cti.full = x_full,
                    qc_full = qc_full)[, ex_ff, drop = FALSE]
    out <- cbind(out, exf)
  }
  if ("coherence" %in% what)
    out$coherence <- .coher(cond, std = TRUE, cti = x, qc_ct = qc_ct, names = FALSE)
  if ("redundant" %in% what){
    out$redundant <- red
  }
  out
}


# Aux function clarify_details()
#   what         required measures
#   measures     measures that are meaningful in the the present context
#                (e.g. "coherence" and "redundant" are not for msc's and asf's
#   notrequired  measures that are not needed at asf level in order to be determined for csf
#   available    available measures (when read from a cna-object)
#   warn         Logical: warn if some required measure is not available
clarify_details <- function(what,
    measures = c("inus", "cyclic", "exhaustiveness", "faithfulness", "coherence", "redundant"),
    available = NULL, notrequired = c("cyclic", "coherence", "redundant"),
    warn = TRUE){
  av_null <- is.null(available)
  available <- union(available, notrequired)
  if (is.logical(what)){
    if (!av_null) measures <- intersect(measures, available)
    return(if (what) measures else character(0))
  }
  what <- measures[pmatch(what, measures, nomatch = 0)]
  if (warn && !av_null && length(not_av <- setdiff(what, available)))
    warning("Measure(s) not available: ", paste0(not_av, collapse = ", "),
            ".\nYou may have to specify details in cna().",
            call. = FALSE)
  if (!av_null) what <- intersect(what, available)
  intersect(measures, what)
}
