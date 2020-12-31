
# build dnf 
#  Note that ct is supposed to be a "full" ct
make.dnf <- function(expr, ct, disj = "+"){
  if (expr %in% c("0", "1")) return(expr)
  selC <- selectCases(expr, ct)
  if (nrow(selC) == 0L) return("0")
  if (nrow(selC) == nrow(ct)) return("1")
  sc <- ctInfo(selC)$scores
  terms <- split(rep(colnames(sc), each = nrow(sc))[sc == 1],
                 row(sc)[sc == 1])
  conj <- C_mconcat(terms, "*")
  conj <- conj[is.minimal(conj)]
  C_concat(conj, disj)
}
# Minimize a single condition
#  Note that x is supposed to be a "full" ct
.minim1 <- function(cond, x, maxstep = c(4, 4, 12)){
  cond <- make.dnf(cond, x)
  if (cond %in% c("0", "1")) return(cond)
  y <- as.vector(qcond_bool(cond, ctInfo(configTable(x))$scores))
  if (isConstant(y)) return(as.character(y[[1]]))
  x$..RESP.. <- y
  suppressMessages({
    .cna <- cna(x, ordering = list("..RESP.."), strict = TRUE, 
                maxstep = maxstep, rm.const.factors = FALSE, rm.dup.factors = FALSE)
  })
  .asf <- asf(.cna)
  if (attr(x, "type") == "mv"){
    .asf <- subset(.asf, .asf$outcome == "..RESP..=1")
  }
  lhs(.asf$condition)
}
# Minimize multiple conditions
minimalize <- function(cond, x = NULL, maxstep = c(4, 4, 12)){
  cond <- noblanks(cond)
  if (is.null(x)){
    x <- full.ct(cond)
  } else {
    x <- full.ct(x)
  }
  cti <- ctInfo(x)
  out <- vector("list", length(cond))
  names(out) <- cond
  if (length(cond) == 0) return(out)
  # check for disjunctive normal form
  dnf <- checkValues(cond, c("+", "*"), colnames(cti$scores))
  cond1 <- cond
  cond1[!dnf] <- vapply(cond[!dnf], make.dnf, x,
                        FUN.VALUE = character(1), USE.NAMES = FALSE)
  out[] <- lapply(cond1, .minim1, x = x, maxstep = maxstep)
  if (any(noOutput <- lengths(out) == 0))
    warning("No minimal solution found for condition(s):\n",
            paste0("  ", cond[noOutput], "\n"),
            "You may try to increase maxstep.", call. = FALSE)
  out
}
