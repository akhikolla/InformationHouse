
is.inus <- function(cond, x = NULL, csf.info = FALSE){
  .inus(x, noblanks(cond), csf.info = csf.info)
}

# ==== .inus() ====
# Switch order of first 2 args to provide dispatching on x
# identifies the asf that are INUS within some 1-sided formulas
# Generic function
.inus <- function(x, cond, ...) UseMethod(".inus")


# aux fun
constantCols <- function(x){
  if (nrow(x)) return(matrixStats::colAlls(x == rep(x[1, ], each = nrow(x))))
  rep(TRUE, ncol(x))
}

# ==== Method for class 'cti' ====
#   x         cti
#   cond      character vector with the csf
#   full      Logical: whether x can be considered to be complete.
#             If not, full.ct() will be applied to x
#   const.ok  If FALSE, the output will be FALSE for cond's evaluating to constant value
# value: A logical vector of same length as cond
# See examples
.inus.cti <- function(x, cond, full = FALSE, const.ok = FALSE, csf.info = FALSE){
  if (!full) x <- full.ct(x)

  qtype <- .qcondType(cond, colnames(x$scores), x$type) 
  # ok <- qtype %in% c("constant", "stdBoolean", "stdAtomic", "stdComplex")
  out <- rep(NA, length(cond))
  names(out) <- cond

  # condition type "constant"
  if (any(selCond <- !is.na(qtype) & qtype == "constant")){
    out[selCond] <- .inusInternal_const(cond[selCond], const.ok = const.ok)
    if (all(!is.na(out) & out)) return(out)
  }

  # condition type "stdBoolean"
  if (any(selCond <- !is.na(qtype) & qtype == "stdBoolean")){
    out[selCond] <- .inusInternal_bool(cond[selCond], x, const.ok = const.ok)
    if (all(!is.na(out) & out)) return(out)
  }

  # condition type "stdAtomic"
  if (any(selCond <- !is.na(qtype) & qtype == "stdAtomic")){
    out[selCond] <- .inusInteral_asf(cond[selCond], x)
    if (all(!is.na(out) & out)) return(out)
  }

  # condition type "stdAtomic"
  if (any(selCond <- !is.na(qtype) & qtype == "stdComplex")){
    out_csf <- .inusInteral_csf(cond[selCond], x, info = csf.info)
    out[selCond] <- out_csf
    attr(out, "csf.info") <- attr(out_csf, "csf.info")
  }

  out
}

# Internal function for conds of type "constant"
.inusInternal_const <- function(cond, const.ok = FALSE){
  out <- rep(TRUE, length(cond))
  if (!const.ok) out[cond %in% c("0", "1")] <- FALSE
  out
}
# Internal function for conds of type "stdBool"
.inusInternal_bool <- function(cond, cti, const.ok = FALSE){    
  out <- rep(TRUE, length(cond))
  
  evalCond0 <- qcond_bool(cond, cti$scores)[, , drop = FALSE]
  if (!const.ok){
    # constant?
    out[constantCols(evalCond0)] <- FALSE
    if (!any(out)) return(out)
  }
  
  notConst <- out
  w <- which(notConst)
  disjs <- strsplit(cond[notConst], "+", fixed = TRUE)
  evalDisj <- qcond_bool(unlist(disjs), cti$scores)
  indList <- C_relist_Int(seq_len(ncol(evalDisj)), lengths(disjs))
  conjs <- C_relist_List(strsplit(unlist(disjs), "*", fixed = TRUE), lengths(disjs))
  
  for (i in seq_along(cond[notConst])){
    evalDisj_i <- evalDisj[, indList[[i]], drop = FALSE]
    conjs_i <- conjs[[i]]

    # Is any conjunction redundant?  -> If yes: condition is not INUS
    if (length(conjs_i) > 1L && any(C_redund(!evalDisj_i))){
      out[w][i] <- FALSE
      next
    }

    # Is a factor in a conjunction redundant? -> If yes: condition is not INUS
    evalCond0_i <- evalCond0[, w[i], drop = TRUE]
    if (any(ll <- (lengths(conjs_i) > 1L))){
      brk <- FALSE
      for (j in seq_along(conjs_i)){
        if (!ll[[j]]) next
        if (any(C_redund(cti$scores[evalCond0_i == 0, conjs_i[[j]], drop = FALSE]))){
          out[w][i] <- FALSE
          brk <- TRUE
        }
        if (brk) break
      }
    }
  }
  
  out
}
# Internal function for conds of type "stdAtomic" (asf)
#  * lhs inus
#  * rhs a factor value
#  * outcome is not part of the lhs
#  * there is >=1 case & < nrow(full.ct) 
.inusInteral_asf <- function(cond, cti){
  out <- rep(NA, length(cond))
  names(out) <- cond
  .rhs <- unname(vapply(cond, rhs, character(1)))
  .lhs <- unname(vapply(cond, lhs, character(1)))
  vals <- colnames(cti$scores)
  valsGrouped <- C_relist_Char(vals, cti$nVal)
  names(valsGrouped) <- names(cti$nVal)
  # 1) check that lhs is inus
  ok1 <- .inusInternal_bool(lhs(cond), cti)
  # 2) check that rhs a factor value
  ok2 <- ifelse(.rhs %in% vals, TRUE, NA)
  if (any(is.na(ok2))){
    message("Expecting a single factor value on rhs of an asf - is.inus() returns NA otherwise.")
  }
  # 3) outcome is not part of the lhs
  lhsVals <- hstrsplit(.lhs, c("+", "*"), fixed = TRUE, relist = FALSE)
  getVar <- function(x) rep(names(valsGrouped), lengths(valsGrouped))[match(x, vals)]
  lhsVars <- hrelist(getVar(lhsVals), 
                     collapse(attr(lhsVals, "lengths"), 2))
  ok3 <- m_all(mapply("!=", getVar(.rhs), lhsVars, SIMPLIFY = FALSE))
  # 4) there is >=1 case & < nrow(full.ct) 
  ok4 <- !apply(qcond_asf(cond, sc = cti$scores, force.bool = TRUE), 2, isConstant)
  out[] <- ok1 & ok2 & ok3 & ok4
  out
}

# Internal function for conds of type "stdComplex" (csf)
#   * all asf inus
#   * there is >=1 case & < nrow(full.ct) 
#   * struct redundancies
#   * partial str red
.inusInteral_csf <- function(cond, cti, info = FALSE){
  # asfs are inus
  asfs <- extract_asf(cond)
  asf_not_inus <- !m_all(lapply(asfs, .inusInteral_asf, cti))
  # structural redundancies?
  redund <- m_any(.redund.cti(cti, cond, simplify = FALSE)) 
  # partial structural redundancies?
  redundInCsf <- partiallyRedundant(asfs, cti)
  # constant Factor?
  constFact <- constFact(cond, cti)
  # constant?
  const <- apply(qcond_csf(cond, cti$scores, force.bool = TRUE), 2, isConstant)
  # multiple outcomes?
  multOutcome <- vapply(happly(asfs, function(x) toupper(rhs(x))), 
                        anyDuplicated, 
                        integer(1)) > 0L
  # mv: All values of some variable as outcome?
  if (cti$type == "mv" && any(!multOutcome)){
    multOutcome[!multOutcome] <- allValuesOutcome(asfs[!multOutcome], cti)
  }
  # Output
  out <- rep(TRUE, length(cond))
  names(out) <- cond
  out[asf_not_inus | redund | redundInCsf | constFact | const | multOutcome] <- FALSE
  out[is.na(asf_not_inus) | is.na(redund) | is.na(redundInCsf) | is.na(constFact) | is.na(const) | is.na(multOutcome)] <- NA
  if (info){
    .info <- data.frame(row.names = cond)
    .info[, c("asf_not_inus", "redund", "redundInCsf", "constFact", "const", "multOutcome")] <- 
      data.frame(asf_not_inus, redund, redundInCsf, constFact, const, multOutcome)
    attr(out, "csf.info") <- .info 
  }
  out
}

# Aux function: partial redundancy
partiallyRedundant <- function(x, cti){
  out <- logical(length(x))
  l1 <- lengths(x) == 1
  out[l1] <- FALSE
  for (i in which(!l1)){
    for (j in seq_along(x[[i]])){
      sel <- as.vector(qcond_csf(C_concat(paste0("(", x[[i]][-j], ")"), "*"), cti$scores, force.bool = TRUE)) == 1L
      .pR <- !.inus.cti(subset(cti, sel), cond = lhs(x[[c(i, j)]]), full = TRUE, const.ok = TRUE)
      if (.pR){
        out[[i]] <- TRUE
        break
      }
      if (out[[i]]) break
    }
  }
  out
}

# csf induces a constant outcome factor
constFact <- function(cond, cti.full){
  if (!length(cond)) return(logical(0))
  cond <- noblanks(cond)
  sc <- cti.full$scores
  qc <- qcond_csf(cond, sc, force.bool = TRUE)
  dimnames(qc) <- NULL
  out <- rep(FALSE, length(cond))
  first <- apply(qc, 2, match, x = 1L)
  x <- sc[, cti.full$resp_nms, drop = FALSE]
  outcomes <- happly(extract_asf(cond), rhs)
  xsplit <- split.default(x, col(x))
  for (.col in split.default(x, col(x))) {
    firstVal <- .col[first]
    foundConst <- !colAnys(qc == 1 & .col != matrix(firstVal, 
        nrow = nrow(qc), ncol = ncol(qc), byrow = TRUE))
    out[foundConst] <- TRUE
  }  
  # Below is a version that tests constancy only for factors that act as an outcome in the respective csf
  # As always testing for all factors is simpler and usually faster, the selection of 
  # outcome factors can be omitted.
  # for (i in seq_along(xsplit)){
  #   .col <- xsplit[[i]]
  #   nm <- colnames(x)[[i]]
  #   sel <- m_any(happly(outcomes, `==`, nm))
  #   if (!any(sel)) next
  #   firstVal <- .col[first[sel]]
  #   foundConst <- !colAnys(qc[, sel, drop = FALSE] == 1 & .col != 
  #                            matrix(firstVal, nrow = nrow(qc), ncol = sum(sel),
  #                                   byrow = TRUE))
  #   out[foundConst] <- TRUE
  # }  
  out  
}


# Aux function: All values of some variable as outcome in some csf?
allValuesOutcome <- function(asfs, cti){
  outcomes <- happly(asfs, rhs, relist = FALSE)
  ii <- rep(seq_along(asfs), lengths(asfs))
  outcomes <- outcomes[order(ii, outcomes)]
  v <- substr(unlist(outcomes), 1, 1)
  qwer <- C_mconcat(C_relist_Char(outcomes, rle(paste(ii, v, sep = "."))$lengths), "|")
  asdf <- C_mconcat(with(cti, C_relist_Char(resp_nms, nVal)), "|")
  nUnique <- function(x) length(unique.default(x))
  m_any(relist1(qwer %in% asdf, vapply(C_relist_Char(v, lengths(asfs)), nUnique, integer(1))))
}



# ==== Method for class 'configTable' ====
#   x         configTable
#   cond      character vector with the csf
# value: A logical vector with same length as cond
.inus.configTable <- function(x, cond, ...){
  cti <- ctInfo(x)
  .inus.cti(cti, cond, ..., full = FALSE) 
}

# ==== Default Method (for matrix or data.frame) ====
.inus.default  <- function(x, cond, ...){
  if (is.null(x)){
    x <- full.ct(cond)
    if (attr(x, "type") == "mv") message("is.inus() with cond of type \"mv\" usually requires explicit specification of x")
  } else {
    x <- full.ct(x)
  }
  cti <- ctInfo(x)
  .inus.cti(cti, cond, ..., full = TRUE)
}
