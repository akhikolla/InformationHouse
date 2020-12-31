#' @rdname values
#' @param ... Arguments passed to class-specific methods. (Not used.)
#' @usage values(x,...)
"values" 

if(requireNamespace("raster", quietly = TRUE)){
	setGeneric("values", def=raster::values)
}else{
	setGeneric(
		name="values",
		def=function(x,...){
			standardGeneric("values")
		}
	)
}

#' Extract and replace values from a gridlayer-derived object (e.g. \code{link{facelayer}}). 
#'
#' The function will get the \code{@values} slot of a \code{\link{facelayer}} object.
#'
#' @param x (\code{\link{facelayer}}) Object.
#' @rdname values
#' @exportMethod values
setMethod(	
	f="values",
	signature="gridlayer",
	definition= function(x){
		return(x@values)
	}
)


#' @usage values(x) <- value
#' @rdname values
"values<-"

if(requireNamespace("raster", quietly = TRUE)){
	setGeneric("values<-", def=raster::`values<-`)
}else{
	setGeneric(
		name="values<-",
		def=function(x,value){
			standardGeneric("values<-")
		}
	)
}
	

#' @param value (\code{logical}, \code{character} or \code{numeric}) Replacement values.
#' @rdname values
#' @exportMethod values<-
setReplaceMethod(	
	f="values",
	signature="gridlayer",
	definition= function(x, value){
		x@values<-value
		if(length(x@values)!=x@length){
			stop("Wrong replacement length.")
		}else{
			return(x)
		}
	}
)


#' The number of faces in an icosahedral grid.
#' 
#' @rdname length
setMethod(	
	f="length",
	signature="gridlayer",
	definition= function(x){
		x@length
	}
)

#' The face names in a \code{\link{facelayer}} class object
#'
#' Function to extract the registered face names to which the \code{\link{facelayer}} renders information.
#'
#' @param x (\code{\link{facelayer}}) Object.
#' @return A vector of \code{character} values, the names of the faces.
#' @rdname names
#' @exportMethod names
setMethod(	
	f="names",
	signature="gridlayer",
	definition= function(x){
		x@names
	}
)


