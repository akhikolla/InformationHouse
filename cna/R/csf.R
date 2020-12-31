
# Extract csf from cna object
csf <- function (x, n.init = 1000, details = x$details,
                 asfx = asf(x, details, warn_details = FALSE),
                 inus.only = x$inus.only, minimalizeCsf = inus.only,
                 acyclic.only = x$acyclic.only, cycle.type = x$cycle.type, 
                 verbose = FALSE){

  # details
  details <- clarify_details(details, available = x$details, 
                             notrequired = c("cyclic", "exhaustiveness", "faithfulness", "coherence", "redundant"), 
                             warn = TRUE)
  details <- union("inus", details)
  if (inus.only && !x$inus.only){
    stop("Calling csf() with inus.only=TRUE on a cna-object generated as cna(..., inus.only = FALSE) is not meaningful.")
  }
  # Checking/enforcing compatibility of inus.only and minimalizeCsf
  if (inus.only && !minimalizeCsf){
    minimalizeCsf <- FALSE
    message("Calling csf() with inus.only=TRUE and minimalizeCsf=FALSE is not meaningful - setting minimalizeCsf=TRUE")
  }
  if (acyclic.only){
    details <- union(details, c("cyclic"))
  }
  if ("cyclic" %in% details) cycle.type <- match.arg(cycle.type)

  # Output if no cond
  if (nrow(asfx) == 0 || n.init <= 0) {
    out <- emptyCondTbl("stdComplex", details = x$details)
    # if (coh) out$coherence <- numeric(0)
    return(out)
  }
  cti <- ctInfo(x$configTable)
  
  # ---- Get csf's as combinations of asf's ---- 
  splitasf <- split(asfx, asfx$outcome)
  n.asf <- vapply(splitasf, nrow, integer(1))
  n.out <- prod(n.asf)
  if (verbose){
    cat("no. of asf", if (inus.only) " (inus.only)", ": ", 
        paste0(names(n.asf), ":", n.asf, collapse = "; "),
        "\nCombine these to ", paste(n.asf, collapse = "*"),
        " = ", n.out, " initial csfs.\n", sep = ""
            )
  }
  splitasf <- lapply(splitasf, 
    function(x) x[head_with_ties(x$complexity, n.init), , drop = FALSE])
  out <- combineAsf(splitasf, n.init)
  initialCsfs <- out$condition
  if (n.init < n.out){
    out <- head(out, n.init)
    warning("Not all csf solutions have been recorded. csf() with a higher value of n.init might find more solutions.")
    if (verbose){
      cat("Selection of the n.init=", n.init, " best csfs: no. of csf ", 
          n.out, "->", n.out <- n.init, "\n", sep = "")
    }
  }

  # ---- Identify structural redundancies from list solutions and reshape result ----
  if (minimalizeCsf){
    out <- minimalizeCsf_and_reshape(out, x$configTable, 
                                     outcomes = names(splitasf),
                                     details = details)
    if (verbose){
      cat("Reduction of csf with structural redundancies (as minimalizeCsf=TRUE): no. of csf: ", 
          n.out, "->", n.out <- nrow(out), "\n", sep = "")
    }
  }
  
  # ---- Calculate required details if not done in the previous step ----
  if (length(details_without_cyclic <- setdiff(details, "cyclic"))){
    if (!x$inus.only) inus_from_asfs <- out$inus
    out[details_without_cyclic] <- 
      .det.cti(cti, out$condition, 
               what = details_without_cyclic, available = details,
               cycle.type = cycle.type, in.csf = TRUE)[details_without_cyclic]
    if (!x$inus.only) out$inus <- out$inus & inus_from_asfs
  }

  # ---- Remove partial structural redundancies (if inus.only=TRUE) ---
  if (inus.only){
    out <- out[out$inus, , drop = FALSE]
    if (verbose){
      cat("Elimination of csf that are constant or have constant factors or partial structural redundancies (as inus.only=TRUE): no. of csf: ",
          n.out, "->", n.out <- nrow(out), "\n", sep = "")
    }
  }
  
  # ---- Remove cyclic solutions if acyclic.only=TRUE ----
  if ("cyclic" %in% details){
    out$cyclic <- cyclic(out$condition, cycle.type = cycle.type, use.names = FALSE)
    if (acyclic.only){
      out <- out[!out$cyclic, , drop = FALSE]
      if (verbose){
        cat("Elimination of cyclic csf (as acyclic.only=TRUE): no. of csf: ", n.out, "->", n.out <- nrow(out),
            "\n", sep = "")
      }
    }
  }
  
  # --- Adjustments in table of csf's, if required --- 
  # if new csfs have been introduced in minimalizeCsf-step:
  if (minimalizeCsf){
    newCsfs <- !(out$condition %in% initialCsfs)
    if (any(newCsfs)){
      # Compute complexity and inus for csf that resulted from removing struct redundancies
      out$complexity[newCsfs] <- getComplexity(out$condition[newCsfs])
      if (any(calc_inus <- is.na(out$inus))){
        asdf <- extract_asf(out$condition[calc_inus])
        out$inus[calc_inus] <- m_all(C_relist_Log(asfx$inus[match(unlist(asdf), asfx$condition)], lengths(asdf)))
      }
      # check if con/cov-thresholds are still attained
      d.eps <- nrow(x$configTable) * .Machine$double.eps
      con <- max(x$con - d.eps, 0)
      cov <- max(x$cov - d.eps, 0)
      ok <- out$consistency >= con & out$coverage >= cov
      if (!all(ok)){
        out <- out[ok, ]
        if (verbose){
          cat("Recalculating con&cov - ", sum(!ok), 
              " csf not reaching the thresholds removed: no. of csf: ", 
              n.out, "->", n.out <- nrow(out), "\n", sep = "")
        }
      }
      # re-order csf
      ord <- with(out, order(complexity, -consistency*coverage))
      out[] <- out[ord, ]
    }
  }

  # --- Remove paranthesis from csf consisting of a simple asf:
  csf_length1 <- !grepl(")*(", out$condition, fixed = TRUE)
  out$condition[csf_length1] <- sub("^\\((.+)\\)$", "\\1", out$condition[csf_length1])
  rownames(out) <- NULL
  as.condTbl(out, condClass = "stdComplex")
}

# Combine asf into csf's
combineAsf <- function(csflist, n){  
  n.asf <- vapply(csflist, nrow, integer(1))
  n.asfCombs <- prod(n.asf)
  l <- length(csflist)
  a <- csflist[[1]][c("complexity", "consistency", "coverage", "inus")]
  a$id1 <- seq_len(nrow(a))
  if (nrow(a) > n){
    ord <- order(a$complexity)
    a <- a[ord, ]
    a <- a[head_with_ties(a$complexity, n), , drop = FALSE]
  }
  if (l >= 2) for (i in seq(2, l)){
    b <- csflist[[i]][c("complexity", "consistency", "coverage", "inus")]
    names(b) <- paste0(names(b), ".",  1)
    b[[paste0("id", i)]] <- seq_len(nrow(b))
    ab <- expand.frames(a, b)
    ab$complexity <- ab$complexity + ab$complexity.1
    con.smaller <- ab$consistency.1 < ab$consistency
    ab$consistency[con.smaller] <- ab$consistency.1[con.smaller]
    cov.smaller <- ab$coverage.1 < ab$coverage
    ab$coverage[cov.smaller] <- ab$coverage.1[cov.smaller]
    ab$inus <- ab$inus & ab$inus.1
    ord <- order(ab$complexity)
    ab <- ab[ord, , drop = FALSE]
    a <- ab[head_with_ties(ab$complexity, n), , drop = FALSE]
    a$complexity.1 <- a$consistency.1 <- a$coverage.1 <- NULL
  }
  a <- a[with(a, order(complexity, -consistency * coverage)), , drop = FALSE]
  selectedRows <- head_with_ties(cbind(a$complexity, 
                                       with(a, -consistency*coverage)), 
                                 n)
  a <- a[selectedRows, , drop = FALSE]
  id <- a[grepl("^id", names(a))]
  allAsfs <- mapply(function(x, i) x$condition[i], csflist, id,
                    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  allAsfs <- C_relist_Char(do.call(rbind, allAsfs), rep(l, nrow(a)))
  data.frame(
    outcome = C_concat(names(csflist), sep = ","),
    condition = C_mconcat(happly(allAsfs, function(x) paste0("(", x, ")")), 
                          "*"),
    consistency = a$consistency, coverage = a$coverage, complexity = a$complexity, 
    inus = a$inus,
    stringsAsFactors = FALSE)
}

# Auxiliary function expand.frames
expand.frames <- function(x, y){
  nx <- nrow(x)
  ny <- nrow(y)
  cbind(x[rep(seq_len(nx), each = ny), , drop = FALSE],
        y[rep(seq_len(ny), nx), , drop = FALSE])
}

# Auxiliary function head_with_ties
head_with_ties <- function(x, n){
  x <- as.matrix(x)
  if (nrow(x) <= n) {
    n <- nrow(x)
  }
  else {
    notDup <- which(!duplicated(x))
    if (all(notDup <= n)) {
      n <- nrow(x)
    }
    else {
      n <- min(notDup[notDup > n]) - 1L
    }
  }
  seq_len(n)
}

# Identify structural redundancies, list solutions and reshape result
minimalizeCsf_and_reshape <- function(out, ct, outcomes, details){
  out1 <- minimalizeCsf(sub(",", "*", out$condition, fixed = TRUE),
                        ct = ct)
  redundancies <- out1$n.asf < length(outcomes)
  out1$n.asf <- out1$redundantParts <- NULL
  nms <- names(out1)
  names(out1)[match("con", nms)] <- "consistency"
  names(out1)[match("cov", nms)] <- "coverage"
  # Add (known) complexity values for unreduced conditions
  unreduced <- !is.na(.pos <- match(out1$condition, out$condition))
  out1$complexity <- rep(NA_integer_, nrow(out1))
  out1$complexity[unreduced] <- out$complexity[.pos[unreduced]]
  out1$inus <- rep(NA, nrow(out1))
  out1$inus[unreduced] <- out$inus[.pos[unreduced]]
  class(out1) <- c("condTbl", "data.frame")
  class(out1$outcome) <- "outcomeString"
  attributes(out1$condition) <- list(class = c("stdComplex", "character"))
  out1$redundant <- rownames(out1) <- NULL
  out1
}
# Aux fun to read complexity of csf (read from character string)
getComplexity <- function(cond){
  if (length(cond) == 0) return(integer(0))
  lhsides <- lapply(extract_asf(cond), lhs)
  ll <- lengths(strsplit(unlist(lhsides), "[\\+\\*]"))
  vapply(C_relist_Int(ll, lengths(lhsides)), sum, integer(1))
}
