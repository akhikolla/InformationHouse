

hutilscpp_rev <- function(x) {
  if (!length(x)) {
    return(x)
  }
  # to avoid slowness on certain altrep vectors
  if (is_altrep(x)) {
    return(seq.int(x[length(x)], x[1]))
  }
  x[length(x):1L]
}

