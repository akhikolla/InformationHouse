
cna <- function (x, type,
    ordering = NULL, strict = FALSE,
    con = 1, cov = 1, con.msc = con, 
    notcols = NULL, rm.const.factors = TRUE, rm.dup.factors = TRUE,  
    maxstep = c(3, 4, 10), inus.only = only.minimal.msc && only.minimal.asf,
    only.minimal.msc = TRUE, only.minimal.asf = TRUE, maxSol = 1e6, 
    suff.only = FALSE, what = if (suff.only) "m" else "ac",
    cutoff = 0.5, border = c("down", "up", "drop"),
    details = FALSE, acyclic.only = FALSE, cycle.type = c("factor", "value"))
{
  # call and type
  cl <- match.call()
  if (!is.null(attr(x, "type"))) type <- attr(x, "type")
  if (missing(type)) type <- "cs"

  # Checking/enforcing compatibility of inus.only with only.minimal.msc and only.minimal.asf
  if (inus.only){
    if (!only.minimal.msc) 
      stop("Calling cna() with inus.only=TRUE and only.minimal.msc=FALSE is not meaningful.")
    if (!only.minimal.asf) 
      stop("Calling cna() with inus.only=TRUE and only.minimal.asf=FALSE is not meaningful.")
  }
    
  # variable names  
  ct <- configTable(x, type = type, rm.const.factors = rm.const.factors, 
                    rm.dup.factors = rm.dup.factors)
  if (any(!grepl("[[:alpha:]]", names(ct)))) 
    stop("All column names must contain a letter.")
  if (anyDuplicated(tolower(names(ct))) != 0L)
    stop("Inadmissable column names: names must be unique when case is ignored!")
  names(ct) <- toupper(names(ct))

  # more formal requirements  
  if (nrow(ct) <= 1 || ncol(ct) <= 1)
    stop("Config table must have at least two rows and two columns.")
  ordering <- check.ordering(ordering, ct)
  # details
  details <- clarify_details(details)
  if (!suff.only) details <- union("inus", details)
  details1 <- setdiff(details, c("coherence", "redundant", "cyclic"))
  if (acyclic.only){
    details <- union(details, c("cyclic"))
  }

  # notcols and ct.out
  if (!is.null(notcols)){
    if (type == "mv") stop("\"notcols\" not applicable if type==\"mv\"")
    if (length(notcols) == 1L && notcols == "all"){
      notcols <- names(ct) 
      if ("ALL" %in% names(ct)) 
        warning("'notcols=\"all\"' is ambiguous if a variable has name \"ALL\". ",
                "notcols is applied to _all_ variables.", call. = FALSE)
    }
    if (!is.character(notcols)) notcols <- names(ct)[notcols]
    notcols <- toupper(notcols)
    if (!all(notcols %in% names(ct))) stop("Wrong specification of 'notcols'")
    ct.out <- ct
    notcols.nrs <- names(ct) %in% notcols
    names(ct.out)[notcols.nrs] <- tolower(names(ct.out)[notcols.nrs])
    ct.out.df <- as.data.frame(ct.out)
    ct.out.df[notcols.nrs] <- lapply(ct.out.df[notcols.nrs], function(x) 1-x)
    attributes(ct.out.df)[c("names", "row.names", "class", "n", "cases", "type")] <- 
      attributes(ct.out)[c("names", "row.names", "class", "n", "cases", "type")]
    ct.out <- ct.out.df
  } else {
    ct.out <- ct
  }

  # Check con and cov values and reduce them slightly to avoid failing 
  # to find conditions due to rounding issues
  if (length(con) != 1 || con < 0 || con > 1 || 
      length(cov) != 1 || cov < 0 || cov > 1){
    stop("Invalid input for 'con' or 'cov'")
  }
  d.eps <- nrow(ct) * .Machine$double.eps
  con <- max(con - d.eps, 0)
  cov <- max(cov - d.eps, 0)

  # maxstep
  stopifnot(length(maxstep) == 3, maxstep > 0)
  maxstep[1:2] <- pmin(maxstep[1:2], maxstep[3])

  # Define config, uniqueValues, resp_nms, factMat, nVal, valueId, scores
  cti <- ctInfo(ct, cutoff = cutoff, border = match.arg(border))
  vnms <- colnames(cti$scores)
  freqs <- attr(ct, "n")
  
  # setup output object
  sol <- vector("list", length(cti$resp_nms))
  names(sol) <- cti$resp_nms
  
  for (zname in cti$resp_nms){
  
    # Identify minimal sufficient conditions
    # --------------------------------------
    # Initialize minSuff, list of minimal sufficient combinations
    .znm <- sub("(.+)=.+", "\\1", zname)
    poteff <- potential.effects(ct, .znm, ordering, strict)
    if (length(poteff) == 0L) next
    
    y <- cti$scores[, zname, drop = TRUE]
    if (zname %in% notcols){
      y[] <- 1-y
      names(sol) <- sub(paste0("^", zname, "$"), tolower(zname), names(sol))
      zname <- tolower(zname)
    }  
    
    minSuff <- findMinSuff(cti, y, poteff, freqs, con.msc, 
                           maxstep, only.minimal.msc)
    if (all(m_is.null(minSuff))) next
    
    msc <- lapply(minSuff, make.msc, outcome = zname, cti = cti, 
                  details = details1, suff.only = suff.only)
    msc <- do.call(rbind, msc)
    nstars <- gregexpr("*", msc$condition, fixed = TRUE)
    msc$complexity <- lengths(nstars) + 1L - (vapply(nstars, "[[", integer(1), 1L) == -1L)
    msc <- msc[order(msc$complexity, -msc$consistency * msc$coverage, msc$condition), , drop = FALSE]
    rownames(msc) <- NULL
    sol[[zname]] <- list(msc = msc)

    if (suff.only) next

    # Find asf's
    # ----------
    .conjList <- lapply(minSuff, "attr<-", "conCov", NULL)
    noMsc <- m_is.null(.conjList)
    .conjList[noMsc] <- lapply(which(noMsc), function(i) matrix(integer(0), 0, i))

    .sol <- findAsf(cti, y, freqs, con, cov, .conjList, maxSol, maxstep, only.minimal.asf)
    stopifnot(length(.sol) == 0 || any(vapply(.sol, is.intList, logical(1))))
      
    if (length(.sol)){
      asf <- make.asf(cti, zname, .sol, inus.only, details1)
      sol[[c(zname, "asf")]] <- asf
    }
  }

  out <- structure(list(), class = "cna")
  out$call <- cl
  out$x <- x
  out$ordering <- ordering
  out$configTable <- ct
  out$configTable_out <- ct.out
  out$solution <- sol
  out$what <- what
  out$details <- details
  out[c("con", "cov", "con.msc", "inus.only", "acyclic.only", "cycle.type")] <- 
    list(con, cov, con.msc, inus.only, acyclic.only, cycle.type)

  return(out)
}


findMinSuff <- function(cti, y, poteff, freqs, con.msc, maxstep, only.minimal.msc){
  nsteps <- min(length(poteff), maxstep[1L])
  minSuff <- vector("list", nsteps)
  if (!only.minimal.msc) minimList <- vector("list", nsteps)
  
  for (.step in seq_along(minSuff)){
    # Select .step-fold conditions occurring in the data
    if (.step == 1L){
      allkfoldConds <- matrix(unlist(cti$config[poteff], use.names = FALSE),
                              ncol = 1)
    } else {
       allkfoldConds <- findAllSubsets(cti$valueId, .step, match(poteff, colnames(cti$valueId)))
    }
    if (!only.minimal.msc){
      minimList[[.step]] <- rep(TRUE, nrow(allkfoldConds))
    }
    # Eliminate combinations that contain a simpler combination that has been selected in an earlier step
    if (.step >= 2L){
      for (k_ in seq_len(.step-1L))
        if (!is.null(minSuff[[k_]])){
          is.minim <- !hasSubsetInM(allkfoldConds, minSuff[[k_]])
          if (only.minimal.msc){
            allkfoldConds <- allkfoldConds[is.minim, , drop = FALSE]
          } else if (!all(is.minim)){
            minimList[[.step]][!is.minim] <- FALSE
          }
      }
      rm(k_)
    }
    if (nrow(allkfoldConds) == 0) break
    stopifnot(anyDuplicated(allkfoldConds) == 0L)
    # Select sufficient conditions and add them to output list
    cons <- conj_conCov(allkfoldConds, cti$scores, y, freqs)
  if (any(isSuff <- cons[1, ] >= con.msc, na.rm = TRUE)){
      isSuff[is.na(isSuff)] <- FALSE
      minSuff[[.step]] <- structure(allkfoldConds[isSuff, , drop = FALSE], 
                                    conCov = cons[, isSuff, drop = FALSE])
      if (!only.minimal.msc) attr(minSuff[[.step]], "minimal") <- minimList[[.step]][isSuff]
    }  
  }
  minSuff
}
make.msc <- function(x, outcome, cti, details, suff.only){
  if (is.null(x)) return(NULL)
  vnms <- colnames(cti$scores)
  out <- data.frame(
    outcome,
    condition = C_mconcat(C_relist_Char(vnms[t(x)], 
                                        rep(ncol(x), nrow(x))), 
                          sep = "*"),
    setNames(data.frame(t(attr(x, "conCov"))), c("consistency", "coverage")),
    stringsAsFactors = FALSE)
  if (!is.null(attr(x, "minimal"))){
    out$minimal <- attr(x, "minimal")
  } else {
    out$minimal <- TRUE
  }
  details <- setdiff(details, "inus")
  if (length(details)){
    out[details] <- .det.cti(cti, paste0(out$condition, "->", out$outcome), 
                             what = details)
  }
  rownames(out) <- NULL
  out
}

findAsf <- function(cti, y, freqs, con, cov, .conjList, maxSol, maxstep, only.minimal.asf){
  .conSc <- lapply(seq_along(.conjList), function(i) C_conjScore(cti$scores, .conjList[[i]]))
  maxstep1 <- maxstep
  maxstep1[[1]] <- min(maxstep1[[1]], length(.conSc))
  .combs <- allStructs(maxstep1)
  
  # initial step  
  .sol <- list()
  i <- 0L
  # search for asf's
  nn <- vapply(.conSc, ncol, integer(1))
  for (i in seq_along(.combs)){
    .c <- .combs[[i]]
    if (prod(nn[.c]) == 0 || any(tabulate(.c, nbins = length(.conSc)) > nn)) next
    .s <- C_find_asf(.c, .conSc[.c], y, freqs, con, cov, maxSol)
    
    if (nrow(.s) > 0){
      .cands <- lapply(seq_along(.c), function(i){
        .part <- .conjList[[.c[i]]][.s[, i] + 1L, , drop = FALSE]
        unname(split.default(.part, row(.part)))
      })
      .newsol <- do.call(mapply, c(list(list), .cands, list(SIMPLIFY = FALSE)))
      stopifnot(is.recIntList(.newsol))
      if (only.minimal.asf && length(.sol)){ 
        .newsol <- .newsol[C_minimal_old(.newsol, .sol, ignore_equals = FALSE)]
      }

      .sol <- c(.sol, .newsol)
    }
  }
  .sol
}
make.asf <- function(cti, zname, .sol, inus.only, details){
  vnms <- colnames(cti$scores)
  .lhs <- hconcat(.sol, c("+", "*"), f = function(i) vnms[i])
  asf <- qcondTbl_asf(paste0(.lhs, "<->", zname), 
                      cti$scores, cti$freq)  
  asf$condition[] <- .lhs
  n_plus_stars <- gregexpr("[\\*\\+]", asf$condition)
  asf$complexity <- 
    lengths(n_plus_stars) + 1L - (vapply(n_plus_stars, "[[", integer(1), 1L) == -1L)
  asf <- asf[order(asf$complexity, -asf$consistency * asf$coverage), , drop = FALSE]
  if (length(details))
    asf <- cbind(asf, 
                 .det.cti(cti, paste0(asf$condition, "<->", asf$outcome), 
                          what = details))
  if (inus.only) asf <- asf[asf$inus, , drop = FALSE]
  rownames(asf) <- NULL
  asf
}

################################################################################

# versions of cna with fixed type
cscna <- function(...){
  cl <- match.call(cna, sys.call())
  stopifnot(is.null(cl$type))
  cl[[1]] <- quote(cna)
  cl$type <- "cs"
  eval.parent(cl)
}
mvcna <- function(...){
  cl <- match.call(cna, sys.call())
  stopifnot(is.null(cl$type))
  cl[[1]] <- quote(cna)
  cl$type <- "mv"
  eval.parent(cl)
}
fscna <- function(...){
  cl <- match.call(cna, sys.call())
  stopifnot(is.null(cl$type))
  cl[[1]] <- quote(cna)
  cl$type <- "fs"
  eval.parent(cl)
}
