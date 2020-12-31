
# Evaluate an 'disjuntion of conjunctions'
#   condstr  string in 'visible' syntax with 'standard boolean' format, i.e. a disjunction of conjunctions
#   sc       'scores' matrix
qcond_bool <- function(condstr, sc){
  spl <- strsplit(condstr, "+", fixed = TRUE)
  conjs <- as.character(unique(unlist(spl, use.names = FALSE, recursive = FALSE)))
  splc <- strsplit(conjs, "*", fixed = TRUE)
  # Check for presence of the factors in the data
  factors <- unique(toupper(unlist(splc, recursive = FALSE, use.names = FALSE)))
  nms <- colnames(sc)
  factor_ok <- !is.na(match(factors, nms))
  if (any(!factor_ok)) stop("Invalid condition specified.", call. = FALSE)
  n <- nrow(sc)
  ms_conjs <- matrix(vapply(splc, function(x) rowMins(sc[, x, drop = FALSE]),
                            numeric(n)),
                     nrow = n, ncol = length(splc))
  out <- matrix(vapply(spl, 
                       function(x) rowMaxs(ms_conjs[, match(x, table = conjs), drop = FALSE]),
                       numeric(n)),
                nrow = n, ncol = length(condstr))
  colnames(out) <- condstr
  out
}

# ===== Handling of 'asf' strings =====

# Evaluate an 'atomic condition' seperated with <-> or ->
#   condstr  string in 'visible' syntax with 'asf' format
#   sc       'scores' matrix
# value: array of dimension nrow(sc) x 2 x length(condstr),
#   with attributes 'condition' and 'response'
qcond_asf <- function(condstr, sc, force.bool = FALSE){
  l <- length(condstr)
  n <- nrow(sc)
  lr <- strsplit(condstr, "<*->")
  out <- qcond_bool(as.character(unlist(lr, recursive = FALSE, use.names = FALSE)), sc)
  if (force.bool){
    dim(out) <- c(n, 2, l)
    equiv <- grepl("<->", condstr, fixed = TRUE)
    out1 <- numeric(n*l)
    dim(out1) <- c(n, l)
    if (any(equiv)) out1[, equiv] <- out[, 1, equiv] == out[, 2, equiv]
    if (any(!equiv)) out1[, !equiv] <- out[, 1, !equiv] <= out[, 2, !equiv]
    colnames(out1) <- condstr
    return(out1)
  }
  response <- colnames(out)[1:l * 2]
  dim(out) <- c(n, 2, l)
  structure(out, 
            condition = condstr, 
            response = response)
}
# quick version of condTbl for 'atomic condition'
# parametes: as above
# Value: condTbl with columns outcome, condition, consistency, coverage
qcondTbl_asf <- function(condstr, sc, freqs){
  cond <- qcond_asf(condstr, sc)
  cond2condTbl(cond, freqs)
}

# Take an expanded 'qcond_asf' array and return the corresponding condTbl
# cond   output from qcondTbl
cond2condTbl <- function(cond, freqs){
  left <- matrix(cond[, 1, , drop = TRUE], dim(cond)[[1]])
  right <- matrix(cond[, 2, , drop = TRUE], dim(cond)[[1]])
  lrmin <- left
  lrmin[right<left] <- right[right<left]
  Sx <- colSums(left * freqs)
  Sy <- colSums(right * freqs)
  Sxy <- colSums(lrmin * freqs)
  if (nrow(cond) == 0){
    consistency <- coverage <- NA_real_
  } else {
    consistency <- Sxy/Sx
    coverage <- Sxy/Sy
  }
  structure(
    data.frame(outcome = attr(cond, "response"),
               condition = structure(attr(cond, "condition"),
                                     class = c("stdAtomic", "character")),
               consistency, coverage,
               row.names = NULL,
               stringsAsFactors = FALSE),
    class = c("condTbl", "data.frame"))
}

# ===== Handling of 'csf' strings =====

# Evaluate an 'complex condition' 
#   condstr  string in 'visible' syntax with 'csf' format
#   sc       'scores' matrix
#   force.bool
#   freqs    if (force.bool & !flat) and freqs is supplied, an attribute conCov is attached
# Value:
# - if flat = TRUE: Returns an array resulting from applying qcond_asf() to the (flattended) 
# vector of all asf present in 'condstr', with additional attribute csflengths
# - if flat = FALSE: Returns a list with a separate object as returned by qcond_asf() for each csf
qcond_csf <- function(condstr, sc, flat = FALSE, force.bool = FALSE,
                      freqs = NULL){
  n <- nrow(sc)
  asfs <- extract_asf(condstr)
  lengths <- lengths(asfs)
  unlasfs <- as.character(unlist(asfs, use.names = FALSE, recursive = FALSE))
  varray <- qcond_asf(unlasfs, sc, force.bool = force.bool)
  if (force.bool){
    relisted <- if (length(varray)){
      relist1(as.vector(varray), lengths*n)
    } else {
      list()
    }
    out <- vapply(
      relisted, 
      function(x){ rowMins(matrix(x, nrow = n)) },
      numeric(n))
    colnames(out) <- condstr
    return(out)
  } else if (!flat & !is.null(freqs)){
    ctbl <- cond2condTbl(varray, freqs)
    asfCons <- C_relist_Num(ctbl$consistency, lengths)
    asfCovs<- C_relist_Num(ctbl$coverage, lengths)
    conCov <- data.frame(
      consistency = vapply(asfCons, min, numeric(1)),
      coverage = vapply(asfCovs, min, numeric(1)))
    conCov$asfCons <- asfCons
    conCov$asfCovs <- asfCovs
    conCov$outcome <- C_mconcat(C_relist_Char(ctbl$outcome, lengths), sep = ",")
  }
  if (flat){
    return(structure(varray, csflengths = lengths))
  }
  gr <- rep(seq_along(lengths), lengths*2*n)
  out <- lapply(split(as.vector(varray), gr),
                function(x){
                  array(x, c(n, 2, length(x)/(2*n)))
                  })
  names(out) <- condstr
  if (exists("conCov", inherits = FALSE))
    attr(out, "conCov") <- conCov
  out
}
# Quick version of condTbl for 'complex condition'
# parametes: as above
# Value: condTbl with columns outcome, condition, consistency, coverage
qcondTbl_csf <- function(condstr, sc, freqs){
  qc <- qcond_csf(condstr, sc, flat = TRUE)
  ctbl <- cond2condTbl(qc, freqs)
  lengths <- attr(qc, "csflengths")
  gr <- rep(seq_along(lengths), lengths)
  data.frame(
    outcome = C_mconcat(split.default(ctbl$outcome, gr), sep = ","),
    condition = condstr,
    con = vapply(split.default(ctbl$consistency, gr), min, numeric(1)),
    cov = vapply(split.default(ctbl$coverage, gr), min, numeric(1)),
    stringsAsFactors = FALSE
  )
}
# Grouped version of condTbl for csf-strings
groupedCondTbl_csf <- function(condstr, sc, freqs){
  qc <- qcond_csf(condstr, sc, flat = TRUE)
  ctbl <- cond2condTbl(qc, freqs)
  lengths <- attr(qc, "csflengths")
  out <- split(ctbl, 
               rep(seq_along(lengths), lengths))
  names(out) <- condstr
  out <- lapply(out, "rownames<-", NULL)
  out
}

