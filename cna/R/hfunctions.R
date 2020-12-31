
# == Contents ==
# hstrsplit
# hconcat
# happly
# relist1
# hlengths / hrelist / hsubset


# ==== hstrsplit ====
# hierarchically apply strsplit the elements of a character vector
# (x can also be a list of character vectors)
hstrsplit <- function(x, split, fixed = TRUE, relist = TRUE, 
                      split.attr = TRUE, lengths.attr = !relist){
  if (xList <- is.list(x)){ 
    u <- unlist(x, recursive = FALSE, use.names = FALSE)
  } else {
    u <- x
  }
  s <- strsplit(u, split[1], fixed = fixed[[1]])
  ll <- lengths(s)
  if (length(split) > 1){
    if (length(fixed) > 1) fixed <- fixed[-1]
    s <- hstrsplit(s, split[-1], fixed = fixed, relist = relist, 
                   split.attr = split.attr)
  }
  if (!relist){
    s <- unlist(s, recursive = TRUE, use.names = FALSE)
    if (lengths.attr){
      attr(s, "lengths") <- c(list(ll), attr(s, "lengths"))
    }
    return(s)
  }
  if (xList){
    s <- C_relist_List(s, lengths(x))
  }
  if (split.attr) attr(s, "split") <- split
  s
}


# ==== hconcat ====
# Hierarchically concatenate a nested list of (usually character) elements
# f is applied to unlist(x) and is expected to be NULL or to return character values
hconcat <- function(x, split = attr(x, "split"), f = NULL){
  if (length(split) > 1){
    u <- unlist(x, recursive = FALSE, use.names = FALSE)
  } else {
    stopifnot(length(split) == 1)
    if (!is.null(f)) x <- happly(x, f)
    return(C_mconcat(x, split))
  }
  r <- hconcat(u, split[-1], f = f)
  C_mconcat(C_relist_Char(r, lengths(x)), split[[1]])
}


# ==== happly ====
# hierachically apply a function
#   the list is expected to be of homogneous structure, i.e. at each level, 
#   all elements must be of the same type
#   f must be a function that returns an atomic vector of the same length as x.
happly <- function(`_x`, f = function(x) x, ..., relist = TRUE){
  if (!is.list(`_x`)){
    return(f(`_x`, ...))
  }
  u <- unlist(`_x`, recursive = FALSE, use.names = FALSE)
  out <- happly(u, f = f, ...)
  if (!relist) return(out)
  relist1(out, lengths(`_x`))
}

# ==== happly2 ====
# Alternative version of happly
# happly2 <- function(`_x`, f = function(x) x, ..., relist = TRUE){
#   u <- unlist(`_x`, recursive = TRUE, use.names = FALSE)
#   r <- f(u, ...)
#   if (!relist) return(r)
#   hrelist(r, hlengths(`_x`))
# }

# Auxiliary fun hlengths: recursive lengths of a nested list
hlengths <- function(x){
  if (!is.list(x[[1]])) return(list(lengths(x)))
  u <- unlist(x, recursive = FALSE, use.names = FALSE)
  c(list(lengths(x)), hlengths(u))
}

# collapse one or several level of a list of integer vectors created with hlenghts()
collapse <- function(x, level){
  stopifnot(is.list(x), level %in% seq_along(x))
  level <- sort.int(level)
  if (length(level) == 1L && level == 1L) return(x[-1])
  out <- x
  l1 <- level[1L]
  out[[l1]] <- NULL
  if (l1>1) out[[l1-1L]] <- vapply(C_relist_Int(x[[l1]], x[[l1-1L]]), sum, 1L)
  if (length(level) == 1L) return (out)
  collapse(out, level[-1L] - 1L)
}

# Auxiliary fun relist1: relist with "method dispatch" on typeof(x)
relist1 <- function(x, g){
  out <- switch(typeof(x),
                integer =   C_relist_Int(x, g),
                double =    C_relist_Num(x, g),
                logical =   C_relist_Log(x, g),
                character = C_relist_Char(x, g),
                list =      C_relist_List(x, g),
                stop("relist1 is unnable to handle data of type \"", typeof(x), "\"")
                )
  out
}

# Auxiliary fun hrelist: recursive relisting
#   struct usually determined by hlengths
# recursive relisting
hrelist <- function(`_x`, struct = attr(`_x`, "lengths"), f = NULL, ...){
  stopifnot(is.list(struct))
  if (!is.null(f) && is.atomic(`_x`))
    `_x` <- f(`_x`, ...)
  l <- length(struct)
  if (l == 0L) return(`_x`)
  r <- relist1(`_x`, struct[[l]])
  if (length(struct) == 1L) return(r)
  hrelist(r, struct[-l])
}

# hsubset:
# subset an unlisted object with a "lengths"-attribute
# subs must be logical, of the same length as x
hsubset <- function(x, subs){
  lsub <- lsubset(attr(x, "lengths"), subs)
  structure(x[lsub], lengths = attr(lsub, "lengths"))
}
lsubset <- function(l, subs, lengths = list()){
  if (length(l) == 0){
    return(structure(subs, lengths = lengths))
  }
  lsubset(l[-1], rep(subs, l[[1]]), c(lengths, list(l[[1]][subs])))
}

