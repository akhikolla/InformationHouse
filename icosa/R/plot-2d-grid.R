#' Plot method for the \code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}} classes
#' 
#'	This function will invoke the \code{plot} method of the \code{\link[sp]{SpatialPolygons}} class.
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}) The object to be plotted.
#' @param projargs (\code{character} or \code{\link[sp:CRS-class]{CRS}}) A projection string for the transformation fo coordinates.
#' @param ... Arguments passed to the \code{plot} function.
#' @rdname plot
#' @return The function has no return value.
#' @exportMethod plot
"plot"

#' Plot method for the trigrid object class
#' @rdname plot
setMethod(
	"plot",
	signature="trigrid",
	definition=function(x,projargs=NULL,...){
		#if no @sp found
		if(suppressWarnings(is.na(x@sp))){
			stop("Slot @sp is empty. Use newsp() to add a 2d respresentation. ")
		}
		
		#transformation is necessary
		if(!is.null(projargs)){
			#	requireNamespace("rgdal")
			if(requireNamespace("rgdal", quietly = TRUE)){
				if(class(projargs)=="CRS"){
					x@sp<-sp::spTransform(x@sp, projargs)
				}
				if(class(projargs)=="character"){
					x@sp<-sp::spTransform(x@sp, sp::CRS(projargs))
				}
			} else{
				stop("The rgdal package is required to appropriately project this object. ")
			}
		
		}
	
		sp::plot(x@sp, ...)
	#	rect(ytop=90, ybot=-90, xleft=-180, xright=180)
	
	}
)

#' Lines method for the \code{trigrid} and \code{hexagrid} classes
#' 
#' This function will invoke the \code{\link[sp:panel]{sp.lines}} method of the \code{\link[sp]{SpatialPolygons}} class.
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}}) Object.
#' @param projargs (\code{character} or \code{\link[sp:CRS-class]{CRS}}) A projection string for the transformation fo coordinates. 
#' @param ... Arguments passed to the \code{\link[sp:panel]{sp.lines}} method.
#' @rdname lines-methods
#' @return The function has no return value.
#' @exportMethod lines
setMethod(
	"lines",
	signature="trigrid",
	definition=function(x,projargs=NULL,...){
		#if no @sp found
		if(suppressWarnings(is.na(x@sp))){
			stop("Slot @sp is empty. Use newsp() to add a 2d respresentation. ")
		}
		
		#transformation is necessary
		if(!is.null(projargs)){
			if(requireNamespace("rgdal", quietly = TRUE)){
				if(class(projargs)=="CRS"){
					x@sp<-sp::spTransform(x@sp, projargs)
				}
				if(class(projargs)=="character"){
					x@sp<-sp::spTransform(x@sp, sp::CRS(projargs))
				}	
			} else{
				stop("The 'rgdal' package is required to appropriately project this object. ")
			}
				
		}
		
		# replace with sp::sp.lines()
		sp::sp.lines(x@sp,...)
	
	}
)

#' Labels of grid vertices, faces and edges.
#' 
#' This function will show where the grid elements are located.
#' @param gridObj (\code{\link{trigrid}}, \code{\link{hexagrid}}) An icosahedral grid.
#' @param type (\code{character}) The type of element to be plotted: either \code{"f"} (faces), \code{"v"} (vertices) or  \code{"e"} (edges).
#' @param projargs (\code{character} or \code{\link[sp:CRS-class]{CRS}}) A projection string for the transformation fo coordinates.
#' @param ... Arguments passed to the \code{\link[graphics]{text}} function.
#' @return The function has no return value.
#' @export
#' @examples
#' gr <- hexagrid(sp=TRUE)
#' plot(gr)
#' gridlabs(gr)
gridlabs<-function(gridObj,type="f",projargs=NULL,...){
	# center back to origin if not there already
		if(gridObj@center[1]!=0 | gridObj@center[2]!=0 | gridObj@center[3]!=0){
			gridObj<-translate(gridObj,-gridObj@center)
		}
		
	if(type=="f"){
		texts<-rownames(gridObj@faceCenters)
		coords<-CarToPol(gridObj@faceCenters, norad=TRUE, origin=gridObj@center)
	}
	
	if(type=="v"){
		texts<-rownames(gridObj@vertices)
		coords<-CarToPol(gridObj@vertices, norad=TRUE, origin=gridObj@center)
	}
	if(type=="e"){
		texts<-rownames(gridObj@edges)
		coord3d<-t(apply(gridObj@edges, 1, function(x){
			apply(gridObj@vertices[x,],2,mean)
		}))
		coords<-CarToPol(coord3d, norad=TRUE, origin=gridObj@center)
		
	}
	
	spPoints<-sp::SpatialPoints(coords, proj4string=CRS("+proj=longlat +a=6371007 +b=6371007"))
	
	#transformation is necessary
	if(!is.null(projargs)){
		#	requireNamespace("rgdal")
		if(requireNamespace("rgdal", quietly = TRUE)){
			if(class(projargs)=="CRS"){
				spPoints<-sp::spTransform(spPoints, projargs)
			}
			if(class(projargs)=="character"){
				spPoints<-sp::spTransform(spPoints, sp::CRS(projargs))
			}

		} else{
			stop("The rgdal package is required to appropriately project this object. ")
		}
		
	}
	
	graphics::text(labels=texts, x=spPoints@coords[,1], y=spPoints@coords[,2],...)
	
}

			


#' Locate grid faces based on their positions on a map
#' 
#' The function returns which grid faces contain the points clicked in a plot.
#' 
#' @param gridObj (\code{\link{trigrid}} or \code{\link{hexagrid}}) The grid object.
#' @param n (\code{integer}) The number of points to be looked up.
#' @param output (\code{character}) Type of output: \code{"faces"} returns only the face names of the points, \code{"full"} returns the coordinates as well.
#' @param ... Arguments passed to the \code{\link[graphics]{locator}} function.
#' 
#' @export
#' @return A vector of \code{character} values, each corresponding to a face identifier.
cellocator <- function(gridObj,n, output="faces",...){
	pointset<- locator(n=n, ...)
	pointset <-cbind(pointset$x, pointset$y)
	cells <- locate(gridObj, pointset)

	if(output=="full"){
		retVal<- data.frame(pointset, cells, stringsAsFactors=FALSE)
	}
	if(output=="faces"){
		retVal<-cells
	}
	return(retVal)
}

