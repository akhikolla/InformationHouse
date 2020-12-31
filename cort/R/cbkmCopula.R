#' @include generics.R empiricalCopula.R
NULL

############################### Checkerboard with known margin copula class #######
.cbkmCopula = setClass(Class = "cbkmCopula", contains = "empiricalCopula",
  slots = c(m = "numeric", margins = "numeric", known_cop = "ANY"),
  validity = function(object) {
   errors <- c()
   if (any(sapply(object@m,function(m){(nrow(object@data)%%m) != 0}))) {
     errors <- c(errors, "m should divide the number of row")
   }
   if (length(object@margins) != dim(object@known_cop)) {
     errors <- c(errors, "number of known margins hsuld be equal to dimention of known copula")
   }
   if(length(object@m) != ncol(object@data)){
     errors <- c(errors, "lengths of m parameter shoudl be the same as the number of columns of data")
   }
   if (length(object@margins) > ncol(object@data)) {
     errors <- c(errors, "the number of known margins should be smaller than the number of total margins in the empirical data !!")
   }
   if (!all(object@margins %in% 1:ncol(object@data))) {
     errors <- c(errors, "provided margins number should be smaller than the number of dimention of empirical data")
   }
   if (length(errors) == 0)
     TRUE else errors
  })
#' cbkmCopula contructor
#'
#' Given some empirical data, and given some known copula estimation on a sub-vector of this data,
#' the checkerboard with known margins construction consist in
#' a conditional pattern where a checkerboard copula is fitted (similar the the `cbCopula` algorithm), but conditionally on some known margins.
#'
#' See the corresponding vignette for more details.
#'
#'
#' @param x the data to be used
#' @param m checkerboard parameter
#' @param pseudo Boolean, defaults to `FALSE`. Set to `TRUE` if you are already providing pseudo-data into the `x` argument.
#' @param margins_numbers numeric integers which determines the margins for the known copula.
#' @param known_cop Copula a copula object representing the known copula for the selected margins.
#'
#' @name cbkmCopula-Class
#' @title Checkerboards with known margins
#' @rdname cbkmCopula-Class
#'
#' @return An instance of the `cbkmCopula` S4 class. The object represent the fitted copula and can be used through several methods to query classical (r/d/p/v)Copula methods, etc.
#' @export
#'
#' @examples
#' dataset <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
#' known_copula <- cbCopula(dataset[,2:3],m=10)
#' (cop <- cbkmCopula(x = dataset,
#'                   m = 5,
#'                   pseudo = TRUE,
#'                   margins_numbers = c(2,3),
#'                   known_cop = known_copula))
cbkmCopula = function(x, m = rep(nrow(x),ncol(x)), pseudo = FALSE, margins_numbers = NULL, known_cop = NULL) {

  # Some error handlers :
  if(missing(x))   { stop("The argument x must be provided") }
  if(ncol(x) == 0) { stop("you are providing a matrix equal to NULL") }
  if(nrow(x) == 0) { stop("you should provide data...") }
  if((is.null(known_cop) && (!is.null(margins_numbers))) || (is.null(known_cop) && (!is.null(margins_numbers)))) { stop("known_cop argument and margins argument must both be provided.") }
  if(!pseudo) { x <- apply(x, 2, rank, na.last = "keep")/(nrow(x) + 1) }
  if(length(m) == 1){m = rep(m,ncol(x))}
  if(length(m) != ncol(x)){stop("You should provide m values same lengths as the number of columns in data.")}
  if(all(is.null(known_cop), is.null(margins_numbers))) { return(.cbCopula(data = as.matrix(x), m = m)) }
  if(!all(margins_numbers %in% 1:ncol(x))){ stop("Margins number should be inside 1:ncol(x), or NULL.") }


  ######## Returning the object :
  .cbkmCopula(data = as.matrix(x),
              dim = ncol(x),
              m = m,
              margins = margins_numbers,
              known_cop = known_cop)
}
setMethod(f = "show",    signature = c(object = "cbkmCopula"),                definition = function(object)    {
  cat("This is a cbkmCopula , with : \n", "  dim =", dim(object), "\n   n =",
      nrow(object@data), "\n   m =", object@m, "\n")
  cat("The variables names are : ", colnames(object@data), "\n")
  cat("The variables ", object@margins, " have a known copula  given by :\n")
  writeLines(paste("\t", capture.output(show(object@known_cop)), sep = ""))
})

#' @describeIn rCopula-methods Method for the cbCopula
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "cbkmCopula"), definition = function(n, copula) {

  # if n=0, return a 0xdim(copula) matrix :
  if (n == 0) {return(matrix(NA, nrow = 0, ncol = dim(copula)))}

  # get copula infos :
  J <- copula@margins
  d = dim(copula)
  p = length(J)
  m = copula@m
  seuil_inf <- boxes_from_points(copula@data,m)

  # First step : simulate the known copula model.
  rez <- matrix(NA,nrow=n,ncol=d)
  rez[,J] <- rCopula(n, copula@known_cop)

  # Second step : Calculate the boxes that corespond to those simulations
  # find out wich boxes were simulated on the J part :
  simu_boxes <- boxes_from_points(rez[, J,drop=FALSE],m[J])
  # get the corresponding boxes number
  inf_seuil <- apply(simu_boxes,1,function(x){
    #
    possibles_boxes = seuil_inf[colSums(t(seuil_inf[,J]) == x) == p,,drop=FALSE]
    if(nrow(possibles_boxes) == 0){
      return(rep(-10,d - p)) # NO boxes were founded -> error code.
    } else {
      return(possibles_boxes[resample(1:nrow(possibles_boxes),size=1),-J])
    }
  })

  # deal with errors :
  errors = colSums(inf_seuil) < 0
  sup_seuil = inf_seuil + 1/m[-J]
  inf_seuil[,errors] <- 0
  sup_seuil[,errors] <- 1

  # simulation :
  rng <- matrix(runif((d - p) * n), nrow = d-p, ncol=n)
  rez[, -J] <- t(inf_seuil + rng * (sup_seuil - inf_seuil))

  return(rez)
})

#' @describeIn pCopula-methods Method for the cbCopula
setMethod(f = "pCopula", signature = c(u = "matrix", copula = "cbkmCopula"),  definition = function(u, copula) {
  # this function implements the formula for the mesure of the copula
  # given in the paper.  remind that pCopula and dCopula generics already
  # transform inputs into matrices...

  ######## Precalculations :
  J               <- copula@margins
  d               <- dim(copula)
  m               <- copula@m
  boxes           <- as.matrix(do.call(function(...){expand.grid(...,KEEP.OUT.ATTRS = FALSE)}, lapply(m,function(x){seq(0,1-1/x,length=x)})))
  colnames(boxes) <- NULL
  rez_frame <- rep(0,nrow(boxes))
  seuil_inf <- boxes_from_points(copula@data,m)
  almost_0  <- 1/(10*max(m))

  ######## Calculations :
  result <- sapply(1:nrow(u),function(n_obs){
    rez          <- rez_frame
    max_box      <- t(pmin(t(boxes) + 1/m,u[n_obs,])) # pmin recycling way forces us to use t()
    ok           <- apply(boxes < max_box,1,all)

    if(sum(ok) == 0){ return(rez) } # dunno what to do in this case.
    mes_lebesgue <- apply((max_box[ok,] - boxes[ok,])[,-J],1,prod)*prod(m[-J])

    nb <- apply(boxes[ok,],1,function(x){
      not_is_zero <- abs(t(seuil_inf)-x)>almost_0
      # < almost_0 => it's zero since they are all multiples of 1/m_i
      # colSums(not_is_zero) == 0 gives a 1 for each point in the box.
      return(c(
        sum(colSums(not_is_zero) == 0), # gives the number of points in the box
        sum(colSums(not_is_zero[J,]) == 0) # gives the number of point in the J part of the box.
      ))
    })

    weights <- nb[1,]/nb[2,]
    weights[nb[2,] == 0] <- 1/prod(m[-J])
    rez[ok]      <- vCopula(boxes[ok,J],max_box[ok,J],copula@known_cop)*weights*mes_lebesgue
    return(rez)
  })
  # return the cumulative sum :
  return(colSums(result))

})




