
ct2df <- function(ct, tt){
  
    # Ensure backward compatibility of argument tt
    if (!missing(tt)){
      warning("Argument tt is deprecated in ct2df(); use ct instead.", 
              call. = FALSE)
      if (missing(ct)) ct <- tt
    }
  
  if (is.null(attr(ct, "n", exact = TRUE))) return(as.data.frame(ct))
  n <- nrow(ct)
  df <- as.data.frame(ct)[rep(seq_len(n), attr(ct, "n")), , drop = FALSE]
  rownames(df) <- unlist(attr(ct, "cases"), use.names = FALSE, recursive = FALSE)
  attributes(df) <- attributes(df)[c("names", "row.names", "class")]
  df
}
