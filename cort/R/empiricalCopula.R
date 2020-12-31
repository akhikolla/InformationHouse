#' @include generics.R
NULL

setClass(Class = "empiricalCopula", contains = c("VIRTUAL"),
   slots = c(data = "matrix",dim="numeric"), validity = function(object) {
     errors <- c()
     if(ncol(object@data) == 0){
       errors <- c(errors, "you are providing a matrix equal to NULL")
     }
     if (prod(apply(object@data, 1:2, is.numeric)) != 1) {
       errors <- c(errors, "the data argument must be a numeric matrix")
     }
     if (prod(object@data <= 1) * prod(object@data >=
                                              0) == 0) {
       errors <- c(errors, "the pseudo-data should be numeric between 0 and 1 (both included)")
     }
     if (nrow(object@data) == 0){
       errors <- c(errors,"data matrix must have at least one row.")
     }
     if (length(errors) == 0)
       TRUE else errors
   })

#' @export
pairs.empiricalCopula <- function(x,...){
   N = nrow(x@data)
   suppressWarnings(graphics::pairs(rbind(rCopula(N,x),x@data),
         lower.panel = function(x,y,ind.lower,col.lower,ind.upper,col.upper,...){suppressWarnings(graphics::points(x[ind.lower],y[ind.lower],col = col.lower,...))},
         upper.panel = function(x,y,ind.upper,col.upper,ind.lower,col.lower,...){suppressWarnings(graphics::points(x[ind.upper],y[ind.upper],col = col.upper,...))},
         ind.upper = c(rep(TRUE,N),rep(FALSE,N)),
         ind.lower = c(rep(FALSE,N),rep(TRUE,N)),
         col.upper = 'red',
         col.lower = 'black',
         pch = 20,
         ...))
}


