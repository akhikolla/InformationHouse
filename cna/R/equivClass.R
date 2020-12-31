
equivClass <- function(conds, x){
  if (missing(x)) x <- conds
  x <- full.ct(x)
  qco <- qcond_bool(conds, ctInfo(x)$scores)
  ec <- apply(qco, 2, C_concat, "|")
  out <- unname(split(conds, ec))
  out <- out[order(lengths(out), decreasing = TRUE)]
  structure(out, class = "equivClass")
}

print.equivClass <- function(x, ...){
  ll <- lengths(x)
  cat("'equivClass'\n",
      sum(ll), " conditions in ", length(ll), " equivalence classes:\n",
      sep = "")
  tbl <- rev(table(ll))
  cs <- unname(cumsum(tbl))
  rg <- as.character(c(1, cs[seq_len(length(cs)-1)]+1))
  rg[tbl>1] <- paste(rg[tbl>1], cs[tbl>1], sep = "-")
  rg <- paste0("(", rg, ")")
  writeLines(paste0("  ", format(tbl), " class", ifelse(tbl>1, "es", ""), 
                    " with ", names(tbl), " elements ", rg))
}

