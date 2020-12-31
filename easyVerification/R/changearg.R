#' Change Function Default Arguments
#' 
#' Override default arguments of functions. This functionality is used to
#' deal with the updated default representation in the SpecsVerification 
#' package (>= v0.5).
#' 
#' @param FUN name of function
#' @param ... arguments to be overriden (e.g. format = 'member')
#' 
#' 
changearg <- function(FUN, ...){
  .FUN <- FUN
  argnames <- names(formals(FUN))
  args <- list(...)
  argi <- which(names(args) %in% argnames)
  if (length(argi) > 0){
    invisible(lapply(argi, function(i) {
      formals(.FUN)[[names(args)[i]]] <<- args[[i]]
    }))
  }
  .FUN
}
