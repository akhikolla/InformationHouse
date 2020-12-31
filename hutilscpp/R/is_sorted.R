


is_sorted <- function(x) {
  stopifnot(is.atomic(x), is.numeric(x))
  switch(typeof(x),
         "integer" = is_sorted_int(x),
         "double" = is_sorted_dbl(x))
}

isntSorted <- function(x) {
  switch(typeof(x),
         "integer" = do_isntSorted_int(x),
         "double" = do_isntSorted_dbl(x))
}
