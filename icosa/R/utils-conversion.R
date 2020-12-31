#' Conversion of polar coordinates to 3d Cartesian coordinates
#' 
#' The function uses basic trigonometric relationships to transform longitude/latitude coordinates on a sphere to xyz Cartesian coordinates.
#' 
#' The authalic mean radius of Earth (6371.007 km) is used by this function as a default. The origin is \code{c(0,0,0)}. The precision of these conversions is not exact (see example \code{c(0,90)} below),
#' but should be considered acceptable when applied at a reasonable scale (e.g. for global analyses using data above \code{10e-6} meters of resolution).
#' 
#' @param x (\code{matrix}, \code{numeric}, \code{data.frame}) A 2-column \code{numeric} matrix with the longitude/latitude data.
#' 
#' @param radius (\code{numeric}) The radius of the sphere. Defaults to the R2 radius of Earth (6371.007km).
#'
#' @param origin (\code{numeric}) Vector with length \code{3}, the XYZ coordinates of the sphere center.
#' @param ... Arguments passed to class-specific methods.
#' 
#' @return An xyz 3-column numeric \code{matrix}, \code{data.frame} or \code{numeric}, depending on the class of \code{x}.
#' 
#' @exportMethod PolToCar
#' @rdname PolToCar
#'
#' @examples 
#' longLat <- rbind(
#'   c(0,0),
#'   #note the precision here!
#'   c(0, 90),
#'   c(-45,12)
#' )
#' # matrix-method
#' xyz <- PolToCar(longLat)
#' # numeric-method
#' xyz2 <- PolToCar(longLat[1,])
#' # data.frame method
#' xyz3 <- PolToCar(as.data.frame(longLat))
setGeneric("PolToCar", function(x,...) standardGeneric("PolToCar"))

# Matrix-method of PolToCar()
#
#' @rdname PolToCar
setMethod(
	"PolToCar",
	signature="matrix",
	function(x, radius=authRadius, origin=c(0,0,0)){
	
		#ignore the NAs
		boolNA<-(is.na(x[,1]) | is.na(x[,2]))
		longLat<-x[!boolNA,, drop=FALSE]
		
		if(length(radius)!=1 | !is.numeric(radius)) stop("Invalid \'radius\' value.")
		
		#essential! check whether the long is first and lat is second
		if(max(abs(longLat[,2]))>90) stop("Latitudinal data should be in the second column of the matrix.")
		
		xVar<-cos(longLat[,2]/180*pi)*cos(longLat[,1]/180*pi)
		yVar<-cos(longLat[,2]/180*pi)*sin(longLat[,1]/180*pi)
		zVar<-sin(longLat[,2]/180*pi)
		
		newMat<-matrix(NA, ncol=3, nrow=nrow(x))
		newMat[!boolNA,]<-cbind(xVar,yVar,zVar)
		colnames(newMat)<-c("x", "y", "z")

		endRes<-newMat*radius
		
		# do the translocation if necessary
		endRes[,1]<-endRes[,1]+origin[1]
		endRes[,2]<-endRes[,2]+origin[2]
		endRes[,3]<-endRes[,3]+origin[3]
		
		return(endRes)

	}

)

# Numeric method for PolToCar
#
#' @rdname PolToCar
setMethod(
	"PolToCar", 
	signature="numeric",
	function(x, radius=authRadius, origin=c(0,0,0)){
		if(length(x)!=2) stop("Please provide either a matrix, or two numeric values.")
		x <- matrix(x, ncol=2, nrow=1) 
		vec <- PolToCar(x, radius=radius, origin=origin)[1,]
		return(vec)
	}
)	


# Data.frame method for PolToCar
#
#' @param long (\code{character}) If \code{x} is a \code{data.frame}, then the column used as longitudes. 
#' @param lat (\code{character}) If \code{x} is a \code{data.frame}, then the column used as latitudes.
#' @rdname PolToCar
setMethod(
	"PolToCar",
	signature="data.frame",
	function(x, radius=authRadius, origin=c(0,0,0), long=NULL, lat=NULL){
		if(is.null(long) & !is.null(lat) |
			!is.null(lat) & is.null(lat)) stop("Please provide both longitude and latitude columns.")
		if(!is.null(long) & !is.null(lat)){
			# reursive data.frame method
			# check whether the coluns are actually there
			if(!any(long==colnames(x))) stop("'long' is not a column of 'x'")
			if(!any(lat==colnames(x))) stop("'lat' is not a column of 'x'")

			# call method recursively
			ret <- PolToCar(x[, c(long, lat)], radius=radius, origin=origin)
		}else{
			# the data frame has to have two columns
			if(ncol(x)!=2) stop("Please provide a two-column data.frame, or 'long' and 'lat' arguments.")
			# use the matrix-method (base case of recursion)
			if(!is.numeric(x[,1]) | !is.numeric(x[,2])) stop("One or two columns of the data.frame is/are not numeric.")
			# base case executes as the matrix-method
			mat <- as.matrix(x)
			ret <- PolToCar(mat, radius=radius, origin=origin)
		}
		return(as.data.frame(ret))
	}
)




#' Conversion of 3d Cartesian coordinates to polar coordinates
#' 
#' The function uses basic trigonometric relationships to transform XYZ coordinates to polar coordinates
#' 
#' @param x (\code{matrix}, \code{data.frame}, \code{numeric}) A 3 column data matrix with XYZ coordinates in Cartesian space.
#' @param origin (\code{numeric}) Vector with length \code{3}, the XYZ coordinates of the sphere center.
#' @param norad (\code{logical}). Toggles whether the rho coordinate (distance from origin) should be omitted from the output.
#' @param ... Arguments passed to class-specific methods.
#' 
#' @return A 3-column or 2-column \code{numeric}, \code{matrix} or \code{data.frame} with longitude, latitude and, if set accordingly, radius data.
#' 
#' @examples
#' # some random points
#' xyz <- rbind(
#'   c(6371, 0,0),
#'   c(0, 6371,0),
#'   c(1000,1000,1000)
#' )
#' 
#' # conversions
#'   CarToPol(xyz)
#' @exportMethod CarToPol
#' @rdname CarToPol
setGeneric("CarToPol", function(x,...) standardGeneric("CarToPol"))

#' Matrix-method of CartoPol()
#' @rdname CarToPol 
setMethod(
	"CarToPol",
	signature="matrix",
	function(x, norad=FALSE, origin=c(0,0,0)){
		#ignore the NAs
		boolNA<-(is.na(x[,1]) | is.na(x[,2]) | is.na(x[,3]))
		x<-x[!boolNA,, drop=FALSE]
		
		# center the coordinates to the center of c(0,0,0)
		x[,1]<-x[,1]-origin[1]
		x[,2]<-x[,2]-origin[2]
		x[,3]<-x[,3]-origin[3]
		
		
		#transform back to spherical coordinates
			xSign<-sign(x[,1])
		
			theta<-atan(x[,2]/x[,1])
			phi<-atan(sqrt(x[,1]^2+x[,2]^2)/x[,3])
		
		#transform spherical coordinates to long/lat
			theta<-theta/pi*180
			phi<-phi/pi*180
		
			#convert to lat-long
			long<-rep(NA, length(theta))
			lat<-rep(NA, length(phi))
			
			lat[phi>=0]<-90-phi[phi>=0]
			lat[phi<0]<--90-phi[phi<0]
			
			long[xSign<0 & theta<=0]<-180+theta[xSign<0 & theta<=0]
			long[xSign<0 & theta>0]<--180+theta[xSign<0 & theta>0]
			long[xSign>=0]<-theta[xSign>=0]
		
			rho<-sqrt(x[,1]^2+x[,2]^2+x[,3]^2)
			
			if(sum(is.nan(long))>0)
			{
				long[is.nan(long)]<-0
			
			}
			#the polar problem:
			lat[lat==90 & x[,3]<0]<- -90
			
			
		if(norad){
			matLongLat<-matrix(NA, ncol=2, nrow=length(boolNA))
			matLongLat[!boolNA,]<-cbind(long, lat)
			rownames(matLongLat)<-rownames(x)
			colnames(matLongLat)<-c("long", "lat")
			return(matLongLat)
			
		}else{	
			matLongLat<-matrix(NA, ncol=3, nrow=length(boolNA))
			matLongLat[!boolNA,]<-cbind(long, lat, rho)
			rownames(matLongLat)<-rownames(x)
			colnames(matLongLat)<-c("long", "lat", "rho")
			return(matLongLat)
			
		}
	}

)

#' Numeric method for CarToPol
#'
#' @rdname CarToPol
setMethod(
	"CarToPol", 
	signature="numeric",
	function(x, norad=FALSE, origin=c(0,0,0)){
		if(length(x)!=3) stop("Please provide either a matrix, or two numeric values.")
		x <- matrix(x, ncol=3, nrow=1) 
		vec <- CarToPol(x, norad=norad, origin=origin)[1,]
		return(vec)
	}
)	


#' Data.frame method for CarToPol
#'
#' @rdname CarToPol
setMethod(
	"CarToPol",
	signature="data.frame",
	function(x, norad=FALSE, origin=c(0,0,0)){
		
		# the data frame has to have two columns
		if(ncol(x)!=3) stop("Please provide a three-column data.frame or matrix.")
		# use the matrix-method (base case of recursion)
		if(!is.numeric(x[,1]) | !is.numeric(x[,2])| !is.numeric(x[,3])) stop("One or more columns of the data.frame is/are not numeric.")
		# base case executes as the matrix-method
		mat <- as.matrix(x)
		ret <- CarToPol(mat, norad=norad, origin=origin)
	
		return(as.data.frame(ret))
	}
)


