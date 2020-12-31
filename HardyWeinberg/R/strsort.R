strsort <-
function(s) {
  paste(sort(unlist(strsplit(s, ""))), collapse = "")  
}
