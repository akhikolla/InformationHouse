#' which of three vectors are the elements (all, any) true?
#' @param x,y,z Logical vectors. Either the same length or length-1
#' @param And Boolean. If \code{TRUE}, only indices where all of x, y, z
#' are TRUE are returned; if \code{FALSE}, any index where x, y, z
#' are TRUE are returned.
#' @param anyNAx,anyNAy,anyNAz Whether or not the inputs have \code{NA}.
#' @export

which3 <- function(x, y, z,
                   And = TRUE,
                   anyNAx = anyNA(x),
                   anyNAy = anyNA(y),
                   anyNAz = anyNA(z)) {
  stopifnot(is.logical(x), is.logical(y), is.logical(z))
  check_TF(And)
  if (anyNAx || anyNAy || anyNAz) {
    cpp_list <- do_which3(x, y, z, And)
    {cpp_list[[2]]}[seq_len(cpp_list[[1]])]
  } else {
    do_which3_mem(x, y, z, And)
  }
}



