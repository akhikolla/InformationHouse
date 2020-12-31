
# fs2cs: 
# ======
# transforms a fs-ct into a sc-ct
# Note: fs2cs is designed to work with a data.frame, not a configTable!
fs2cs <- function(x, cutoff = 0.5, border = c("down", "up", "drop")){
  if (cutoff <= 0 | cutoff >= 1)
    stop("cutoff must be >0 and <1")
  border <- match.arg(border)
  x.r <- x
  x.r[] <- round(x.r + (0.5 - cutoff))
  if (border == "up") 
    x.r[x%%1 == cutoff] <- x.r[x%%1 == cutoff] + 1
  if (border == "drop") 
    x.r <- x.r[!rowAnys(x%%1 == cutoff), , drop = FALSE]
  x.r
}

# ctInfo: 
# ========
# Auxiliary function that extracts useful informations from a configTable
ctInfo <- function(ct, cutoff = 0.5, border = c("down", "up", "drop")){
  stopifnot(inherits(ct, "configTable"))
  type <- attr(ct, "type")
  border <- match.arg(border)
  
  if (type == "cs") {
    uniqueValues <- lapply(ct, function(x) 1:0)
    resp_nms <- names(ct)
    valueId <- 2 - as.matrix(ct)
  } else if (type == "mv") {
    uniqueValues <- lapply(ct, function(x) sort(unique.default(x)))
    resp_nms <- mapply(paste, names(ct), uniqueValues, MoreArgs = list(sep = "="), 
                       SIMPLIFY = FALSE)
    resp_nms <- unlist(resp_nms, use.names = FALSE)
    valueId <- mapply(match, ct, uniqueValues, SIMPLIFY = TRUE, USE.NAMES = TRUE)
    if (!is.matrix(valueId))
      valueId <- matrix(valueId, nrow = nrow(ct), 
                        dimnames = list(NULL, names(uniqueValues)))
  } else if (type == "fs") {
    ct.r <- fs2cs(as.data.frame(ct), cutoff = cutoff, border = border)
    uniqueValues <- lapply(ct.r, function(x) 1:0)
    resp_nms <- names(ct)
    valueId <- 2 - as.matrix(ct.r)
  }
  stopifnot(resp_nms == toupper(resp_nms))
  nVal <- lengths(uniqueValues)
  valueId <- sweep(valueId, 
                   2, 
                   c(if (length(nVal)) 0, cumsum(nVal[-length(nVal)])), 
                   "+")
  stopifnot(valueId%%1 == 0)
  mode(valueId) <- "integer"
  config <- split(seq_len(sum(nVal)), 
                  rep(seq_along(nVal), nVal))
  names(config) <- names(ct)
  if (prod(dim(ct)) == 0){
    scores <- matrix(numeric(0), 0, 0)
  } else {
    scores <- do.call(cbind, lapply(resp_nms, factMat(type), ct = ct))
  }
  freq <- attr(ct, "n")
  structure(mget(c("type", "resp_nms", "nVal", "uniqueValues", "config", 
                   "valueId", "scores", "freq")),
            class = c("cti", "tti"))
}
# function to build 'scores' matrix
factMat <- function(type){
  switch(type, 
    cs =  function(x, ct = tt, tt){
      structure(cbind(ct[[x]], 1 - ct[[x]]), 
                .Dimnames = list(NULL, c(x, tolower(x))))},
    mv = function(x, ct = tt, tt){
      v <- eval(parse(text = sub("=", "==", x, fixed = TRUE), keep.source = FALSE), ct)
      matrix(as.integer(v), dimnames = list(NULL, x))},
    fs = factMat <- function(x, ct = tt, tt){
      structure(cbind(ct[[x]], 1 - ct[[x]]), 
                .Dimnames = list(NULL, c(x, tolower(x))))}
  )
}
  
# subset method for class "cti"
subset.cti <- function(x, i, ...){
  x$valueId <- x$valueId[i, , drop = FALSE]
  x$scores <- x$scores[i, , drop = FALSE]
  x$freq <- x$freq[i]
  x
}

# old function name (for compatibility of material outside the package)
tt.info <- ctInfo


# function check.ordering
# =======================
check.ordering <- function(ordering, ct){
  if (is.null(ordering)) return(NULL)
  if (!is.list(ordering))
    stop("ordering must be NULL or a list.")
  # set all names to uppercase
  ordering <- happly(ordering, toupper)
  ord.vars <- unlist(ordering, use.names = FALSE)
  if(!all(ord.vars %in% colnames(ct))){
    stop("Factor listed in ordering does not exist: ",
         setdiff(ord.vars, colnames(ct)))
  }
  if (anyDuplicated(unlist(ordering)))
    stop("At least one factors appears twice in the ordering.")
  if (!all(colnames(ct) %in% ord.vars)){
    ordering <- c(list(setdiff(colnames(ct), ord.vars)), ordering)
  }
  ordering
}


# potential.effects:
# ==================
potential.effects <- function (x, zname, ordering = NULL, strict = FALSE){
  poteff0 <- setdiff(names(x), zname)
  if (is.null(ordering)){
    out <- poteff0
  } else {
    rownames(x) <- NULL
    ordLevelZ <- which(vapply(ordering, function(x) any(zname %in% x), logical(1)))
    poteff <- unlist(ordering[seq_len(ordLevelZ - strict)], use.names = FALSE)
    out <- poteff[match(poteff0, poteff, 0)]  # same order as columns in x
  }
  if (attr(x, "type") == "mv"){ # Case mv with expanded configTable... [nonstandard!!]
    vname <- sub("(.+)=.+", "\\1", zname)
    out <- setdiff(out, vname)
  }
  out
}



# intList, intList2:
# ------------------
# An "intList" is a representation of a formula as a list of integers,
# where the integers refer to colnames in scores
# Example: If the colnames in scores are c("A", "a", "B", "b", "C", "c") then
#   "A+b*C" ist then represented as list(1L, c(4L, 5L))
# An "recIntList" is a list of "intList" objects.
is.intList <- function(x){
  is.list(x) && all(vapply(x, is.integer, logical(1)))
}
is.recIntList <- function(x) is.list(x) && all(vapply(x, is.intList, integer(1)))

# All possible structures of asf 
# maxstep integer of length 3:
#   element 1: maximum number of terms in conjunctions
#   element 2: maximum number of conjunctions in disjunction
#   element 3: maximum total number of terms
# examples: 
# struct(c(2, 3, 5))
# struct(c(3, 2, 6))
allStructs <- function(maxstep){
  maxCompl <- maxstep[[3]]
  if (maxCompl == 1) return(list(1L))
  maxnum <- maxstep[[1]]
  maxlen <- maxstep[[2]]
  one.less <- allStructs(c(maxnum, maxlen, maxCompl - 1L))
  out <- c(list(list(rep(1L, maxCompl))),
           lapply(one.less, add12each))
  out <- unlist(out, recursive = FALSE)
  out <- out[vapply(out, max, integer(1)) <= maxnum]
  out <- out[lengths(out) <= maxlen]
  out <- lapply(out, sort)
  unique.default(c(one.less, out))
}
# auxiliary function
add12each <- function(x)
  mapply(function(l, i){
    l[i] <- l[i] + 1L
    l
  },
  rep(list(x), length(x)),
  seq_along(x),
  SIMPLIFY = FALSE)



