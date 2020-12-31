
# rreduce: reduce a cond wrt a data.fame/configTable:
rreduce <- function(cond, x = full.ct(cond), full = !missing(x), 
                    verbose = FALSE, maxiter = 1000,
                    simplify2constant = TRUE){
  if (!inherits(x, "cti")){
    if (!inherits(x, "configTable")){
      x <- configTable(x, rm.dup.factors=FALSE, rm.const.factors=FALSE)
    }
    x <- ctInfo(x)
  }
  stopifnot(length(cond) == 1L, nrow(x) > 1L)
  cond <- noblanks(cond)
  if (x$type == "fs") stop("Invalid use of data of type 'fs'." )
  if (full) x <- full.ct(x)
  sc <- x$scores

  evalCond0 <- drop(qcond_bool(cond, sc))
  if (simplify2constant && all(evalCond0 == evalCond0[1L])) 
    return(as.character(evalCond0[1L]))
  
  disjs <- strsplit(cond, "+", fixed = TRUE)[[1L]]
  evalDisj <- qcond_bool(disjs, sc)
  conjs <- strsplit(disjs, "*", fixed = TRUE)

  counter <- 0L
  
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    # remove redundant disjuncts...
    repeat{
      if (length(conjs) > 1L && any(redundantDisj <- C_redund(!evalDisj))){
        rmDisj <- which(redundantDisj)[sample(sum(redundantDisj), 1L)]
        conjs <- conjs[-rmDisj]
        disjs <- disjs[-rmDisj]
        evalDisj <- evalDisj[, -rmDisj, drop = FALSE]
        changed <- TRUE
        if (verbose) cat(counter, ":", C_concat(disjs, sep = " + "), "\n")
        counter <- counter + 1L
      } else break
    }
    # remove redundant factors...
    repeat {
      if (length(wll <- which(lengths(conjs) > 1L))){
        redundantFactors <- lapply(wll, 
          function(j) C_redund(sc[evalCond0 == 0, conjs[[j]], drop = FALSE]))
        if (any(haveRed <- vapply(redundantFactors, any, logical(1)))){
          pickDisj <- which(haveRed)[sample(sum(haveRed), 1L)]
          redFact <- redundantFactors[[pickDisj]]
          rmFact <- which(redFact)[sample(sum(redFact), 1L)]
          conjs[[wll[[pickDisj]]]] <- conjs[[wll[[pickDisj]]]][-rmFact]
          changed <- TRUE
          if (any(dupDisj <- duplicated(lapply(conjs, sort)))) conjs <- conjs[!dupDisj]
          disjs <- C_mconcat(conjs, sep = "*")
          evalDisj <- qcond_bool(disjs, sc)
          counter <- counter + 1L
          if (verbose) cat(counter, ":", C_charList2string(conjs, disj = " + "), "\n")
        } else break
      } else break
      if (counter > maxiter) stop("Reached maximal number of reduction steps (maxiter)")
    }
  }
  out <- C_charList2string(conjs)
  out
}


# getCond: Define a condition from a data.frame/matrix:
#   "Internal" facility: If asf is NULL, aresult is returned as an 'intList'
#   [used in cnaOpt::cnaOpt]
getCond <- function(x, outcome = NULL, type, asf = TRUE){
  if (!is.null(attr(x, "type"))) type <- attr(x, "type")
  if (missing(type)) type <- "cs"
  if (type == "fs") stop("getCond is not applicable to fs data.")
  ct <- configTable(x, type = type, rm.dup.factors = FALSE, rm.const.factors = FALSE, 
                    verbose = FALSE)
  if (!is.null(outcome)){
    outcomeName <- sub("=.+", "", outcome)
    stopifnot(length(outcome) == 1, outcomeName %in% names(x))
    outcomePositive <- condition(outcome, ct)[[1]] == 1
    xx <- as.matrix(ct)[outcomePositive, setdiff(colnames(x), outcomeName), drop = FALSE]
    xx <- unique(xx)
  } else {
    xx <- as.matrix(ct)
  }
  if (type == "cs"){
    b <- matrix(colnames(xx), nrow(xx), ncol(xx), byrow = TRUE)
    b[xx==0] <- tolower(b[xx==0])
  } else {
    b <- matrix(colnames(xx), nrow(xx), ncol(xx), byrow = TRUE)
    b[] <- paste0(b, "=", xx)
  }
  sol <- unname(split(b, row(b)))
  if (is.null(asf)) return(sol)
  sol <- C_mconcat(sol, "*")
  if (!asf) return(sol)
  C_concat(sol, "+")
}
